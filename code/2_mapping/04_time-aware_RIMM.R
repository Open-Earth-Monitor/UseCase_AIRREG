library(dplyr)
library(sf)
library(stars)
library(tictoc)
library(gstat)
library(terra)
library(pbmcapply)
#setwd("AIRREG/spat_stat_cube")
source("../code/airreg_functions.R")  # @Edzer: only change this path

method = "ta_RIMM"
alpha = 0.9  # smoothing factor

# load data --------------------------------------------------------------------
fa = st_read("focus_area.gpkg")
lcs = bind_rows(readRDS("aq_data_hourly_sensor_community.rds"),
                readRDS("aq_data_hourly_opensensemap.rds")) |> 
  #ilter(neighbors_within_100m == 0) |> 
  group_by(sensor_id) |> 
  arrange(time) |> 
  mutate(
    pm25_smoothed = coalesce(
      zoo::rollapply(pm25, width = 6, FUN = function(x) {
      Reduce(function(a, b) alpha * b + (1 - alpha) * a, x, right = TRUE)
    }, fill = NA, align = "right"), pm25),
    pm25_lag1 = coalesce(lag(pm25, 1), pm25),
    pm25_lag24 = coalesce(lag(pm25, 24), pm25),
    pm25_roll_8h = coalesce(zoo::rollmeanr(pm25, 8, fill = NA), pm25)
  ) |> 
  ungroup() |> 
  select(time, longitude, latitude, starts_with("pm25"))

# saveRDS(lcs, "LowCostSensor_pm25_hourly_Feb23.rds")

static = readRDS("static_predictors.rds") |> 
  select(-Land_Cover)
cams = readRDS("cams_pm25_hourly_stars.rds") |> 
  mutate(cams_pm25 = log(cams_pm25))
era_rh = readRDS("era5_rh_hourly_stars.rds")
era_t2m = readRDS("era5_t2m_hourly_stars.rds")
era_wd = readRDS("era5_wd_hourly_stars.rds")
era_ws = readRDS("era5_ws_hourly_stars.rds")
st_dimensions(era_ws)[["x"]]$point = NULL
st_dimensions(era_ws)[["y"]]$point = NULL
st_dimensions(era_ws)[["time"]]$point = NULL

eea = readRDS("aq_data_hourly_eea.rds") |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) |> 
  st_transform(st_crs(static))

# set up timeseries; collect model metrics -------------------------------------

timeseries = st_get_dimension_values(era_rh, "time")

mfile = paste0("metrics_", method, ".rds")
if (!file.exists(mfile)){
  metrics = list(time=NULL, cv=NULL, RMSE=NULL, R2=NULL, MAE=NULL, MBE=NULL)
} else {
  metrics = readRDS(mfile)
  # diff timeseries & metrics! after reading
  timeseries = timeseries[!timeseries %in% metrics$time]
}

# initial state map (= mean pm2.5)
# static$pm25_state = mean(lcs$pm25_smoothed, na.rm=T)

for (i in 1:length(timeseries)){
  
  # read lcs data for current time step
  t = timeseries[i]
  message(t)
  df_now_sf = get_data_now(t, lcs) |> filter(pm25 > 0)
  
  if (nrow(df_now_sf) < 5){
    message("not enough observations. abort!")
    next
  }

  # combine static and dynamic covariates
  cov = static
  cov$cams_pm25 = cams[,,,i, drop=T]$cams_pm25
  cov$rh = era_rh[,,,i, drop=T]$rh
  cov$t2m = era_t2m[,,,i, drop=T]$t2m
  cov$wd = era_wd[,,,i, drop=T]$wd
  cov$ws = era_ws[,,,i, drop=T]$ws 
  
  vars = c("pm25_smoothed", 
           "pm25_lag1", "pm25_lag24", "pm25_roll_8h")
  dynamic = pbmclapply(vars, nn_interpolation, 
                       sf = df_now_sf, 
                       template=static, 
                       mc.cores = length(vars))
  cov = c(cov, do.call(c, dynamic))
  #plot(merge(cov))
  
  # modeling table with target and covariates
  cov_ext = st_extract(cov, df_now_sf) |> st_drop_geometry()
  df_now_sf = bind_cols(df_now_sf |> select(time, pm25), 
                        cov_ext) |> 
    tidyr::drop_na()
  
  # linear model
  LM = lm(log(pm25) ~ Elevation + Population + Highway + Major_Roads + 
            Minor_Roads + CLC_Nat + rh + t2m + wd + ws + cams_pm25 +
            pm25_lag1 + pm25_lag24 + pm25_roll_8h + pm25_smoothed, data = df_now_sf)
  #summary(LM)
  
  # predict + get residuals
  cov = cov[fa[1,]]
  cov$lm_pred = exp(predict(LM, newdata = cov))
  
  res = exp(LM$residuals)
  res = res[res < 100]
  df_now_sf$res = NA
  df_now_sf$res[as.numeric(names(res))] = res
  
  # residual kriging
  cov$res = ord_krig(x=df_now_sf |> filter(!is.na(res)), 
                     var="res", template = cov, 
                     nmin = 10, nmax = 20)
  
  # result
  cov$pred = cov$lm_pred + cov$res
  write_stars(cov["pred"], paste0("maps/map_", method, "_", as.integer(t), ".tif"))
  
  #plot_pm25_pred(cov)

  # validation metrics
  val = eea |> filter(time == as.integer(t))
  val$pred = st_extract(cov["pred"], val) |> pull("pred")
  
  metrics = bind_rows(metrics, 
                      list(time = as.character(t), cv = "val", 
                           model_measures(val$raw_pm2_5, 
                                          val$pred)) |> purrr::flatten(),
                      list(time = as.character(t), cv = "model", 
                           model_measures(LM$model$`log(pm25)`, 
                                          LM$fitted.values)) |> purrr::flatten()
  )

  print(metrics[(nrow(metrics)-1):nrow(metrics),])
  saveRDS(metrics, mfile)
  message("------------------------------------")
  gc()
}

list.files("maps/", method) |> length()


library(ggplot2)
ggplot(metrics, aes(x=time,color = cv)) +
  geom_point(aes(y=RMSE)) +
  geom_point(aes(y=R2), pch = 3)



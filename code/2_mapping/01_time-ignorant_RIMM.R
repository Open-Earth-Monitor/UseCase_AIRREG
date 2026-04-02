library(dplyr)
library(sf)
library(stars)
library(tictoc)
source("AIRREG/code/airreg_functions.R")

# load data --------------------------------------------------------------------

lcs = bind_rows(readRDS("AIRREG/spat_stat_cube/aq_data_hourly_sensor_community.rds"),
                readRDS("AIRREG/spat_stat_cube/aq_data_hourly_opensensemap.rds")) |> 
  filter(neighbors_within_100m == 0) |> 
  select(time, pm25, longitude, latitude)


static = readRDS("AIRREG/spat_stat_cube/static_predictors.rds") |> 
  select(-Land_Cover)
cams = readRDS("AIRREG/spat_stat_cube/cams_pm25_hourly_stars.rds") |> 
  mutate(cams_pm25 = log(cams_pm25))
era_rh = readRDS("AIRREG/spat_stat_cube/era5_rh_hourly_stars.rds")
era_t2m = readRDS("AIRREG/spat_stat_cube/era5_t2m_hourly_stars.rds")
era_wd = readRDS("AIRREG/spat_stat_cube/era5_wd_hourly_stars.rds")
era_ws = readRDS("AIRREG/spat_stat_cube/era5_ws_hourly_stars.rds")
st_dimensions(era_ws)[["x"]]$point = NULL
st_dimensions(era_ws)[["y"]]$point = NULL
st_dimensions(era_ws)[["time"]]$point = NULL

eea = readRDS("AIRREG/spat_stat_cube/aq_data_hourly_eea.rds") |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) |> 
  st_transform(st_crs(static))

# set up timeseries; collect model metrics -------------------------------------

timeseries = st_get_dimension_values(era_rh, "time")
#metrics = list(time=NULL, cv=NULL, RMSE=NULL, R2=NULL, MAE=NULL, MBE=NULL)
metrics = readRDS("AIRREG/spat_stat_cube/metrics_RIMM.rds")

# diff timeseries & metrics! after reading
timeseries = timeseries[!timeseries %in% metrics$time]

for (i in 1:length(timeseries)){
  
  # read lcs data for current time step
  t = timeseries[i]
  message(t)
  df_now_sf = get_data_now(t, lcs) |> filter(pm25 > 0)
  
  # combine static and dynamic covariates
  cov = static
  cov$cams_pm25 = cams[,,,i, drop=T]$cams_pm25
  cov$rh = era_rh[,,,i, drop=T]$rh
  cov$t2m = era_t2m[,,,i, drop=T]$t2m
  cov$wd = era_wd[,,,i, drop=T]$wd
  cov$ws = era_ws[,,,i, drop=T]$ws 
  
  # modeling table with target and covariates
  cov_ext = st_extract(cov, df_now_sf) |> st_drop_geometry()
  df_now_sf = bind_cols(df_now_sf, cov_ext) 
  
  # linear model
  LM = lm(log(pm25) ~ Elevation + Population + Highway + Major_Roads + Minor_Roads + CLC_Nat + 
            rh + t2m + wd + ws + cams_pm25, data = df_now_sf)
  
  # predict + get residuals
  cov$lm_pred = exp(predict(LM, newdata = cov))
  
  res = exp(LM$residuals)
  df_now_sf$res = NA
  df_now_sf$res[as.numeric(names(res))] = res
  
  # residual kriging
  cov$res = ord_krig(x=df_now_sf |> filter(!is.na(res)), 
                     var="res", template = cov, 
                     nmin = 5, nmax = 20)
  
  # result
  cov$pred = cov$lm_pred + cov$res
  write_stars(cov["pred"], paste0("AIRREG/spat_stat_cube/maps/map_RIMM_", as.integer(t), ".tif"))
  
  #plot_pm25_pred(cov)
  
  # validation metrics
  val = eea |> filter(time == as.integer(t))
  val$pred = st_extract(cov["pred"], val) |> pull("pred")
  
  metrics = bind_rows(metrics, 
                      list(time = as.character(t), cv = "val", 
                           model_measures(val$raw_pm2_5, val$pred)) |> purrr::flatten(),
                      list(time = as.character(t), cv = "model", 
                           model_measures(LM$model$`log(pm25)`, LM$fitted.values)) |> purrr::flatten()
  )
  
  print(metrics[(nrow(metrics)-1):nrow(metrics),])
  saveRDS(metrics, "AIRREG/spat_stat_cube/metrics_RIMM.rds")
  message("------------------------------------")
  gc()
}

list.files("AIRREG/spat_stat_cube/maps/") |> length()


library(ggplot2)
ggplot(metrics, aes(x=time,color = cv)) +
  geom_point(aes(y=RMSE)) +
  geom_point(aes(y=R2), pch = 3)

m2 = metrics
m2$time = as.integer(as.POSIXct(m2$time))

lcs_s = group_by(lcs, time) |> 
  summarise(pm25 = mean(pm25, na.rm=T)) |> 
  inner_join(m2 )

ggplot(lcs_s, aes(x=time)) +
  geom_point(aes(y=pm25), color="blue") +
  geom_point(aes(y=RMSE), color="red") +
  geom_point(aes(y=R2), color="orange", pch = 3)

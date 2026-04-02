library(dplyr)
library(sf)
library(stars)
library(tictoc)
library(gstat)
library(terra)
library(pbmcapply)
#setwd("AIRREG/spat_stat_cube")
source("code/airreg_functions.R")  # @Edzer: only change this path

method = "RF_state_smooth_test_noLog_240"
map_dir = file.path("AQ_maps", method, "tifs")
if (!dir.exists(map_dir)) dir.create(map_dir, recursive = T)

# load data --------------------------------------------------------------------
fa = st_read("AQ_data/focus_area.gpkg")

static = readRDS("supplementary/static_predictors.rds") |> 
  select(-Land_Cover)

# hourly time series for 2023
# t = seq(as.POSIXct("2023-01-01T00:00:00", tz = "UTC", format= "%Y-%m-%dT%H"), 
#         as.POSIXct("2024-01-01T00:00:00", tz = "UTC", format= "%Y-%m-%dT%H"), by="hour") |> as.numeric()
# t = t[-length(t)]


# cams = readRDS("cams_pm25_hourly_stars.rds") |> 
#   mutate(cams_pm25 = log(cams_pm25))
# era_rh = readRDS("era5_rh_hourly_stars.rds")
# era_t2m = readRDS("era5_t2m_hourly_stars.rds")
# era_wd = readRDS("era5_wd_hourly_stars.rds")
# era_ws = readRDS("era5_ws_hourly_stars.rds")
# st_dimensions(era_ws)[["x"]]$point = NULL
# st_dimensions(era_ws)[["y"]]$point = NULL
# st_dimensions(era_ws)[["time"]]$point = NULL

dynamic = readRDS("supplementary/test240_cov.rds")
st_dimensions(dynamic)[["x"]]$point = NULL
st_dimensions(dynamic)[["y"]]$point = NULL
st_dimensions(dynamic)[["time"]]$point = NULL

# set up timeseries
timeseries = st_get_dimension_values(dynamic, "time")

# dynamic = read_stars(list.files("supplementary", "resampled", full.names = T)) |> 
#   setNames(c("cams_pm2p5","rh","t2m","wd","ws")) |>
#   st_set_dimensions(names = c("x", "y", "time"), 
#                     which = 3, values = seq(as.POSIXct("2023-01-01"), by="hour", length.out=8760))
# 
# d_list = list.files("supplementary", "resampled", full.names = T)
# 
# read_dynamic = function(i, d_list){
#   
#   c(rast(d_list[1], lyrs=i), 
#     rast(d_list[2], lyrs=i), 
#     rast(d_list[3], lyrs=i),
#     rast(d_list[4], lyrs=i), 
#     rast(d_list[5], lyrs=i)) |> 
#     st_as_stars() |> 
#     setNames(c("cams_pm2p5","rh","t2m","wd","ws")) 
# }


# set up timeseries

#timeseries = seq(as.POSIXct("2023-01-01"), by="hour", length.out=8760)

eea = readRDS("AQ_data/aq_data_cube_hourly_2023_eea.rds") |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) |> 
  st_transform(st_crs(static))

lcs = readRDS("AQ_data/lcs_analysis_ready_data_pm25_hourly_2023.rds") |> 
  filter(time %in% as.integer(timeseries))

# collect model metrics -------------------------------------

mfile = paste0("AQ_maps/metrics_", method, ".rds")
if (!file.exists(mfile)){
  metrics = list(time=NULL, cv=NULL, RMSE=NULL, R2=NULL, MAE=NULL, MBE=NULL)
} else {
  metrics = readRDS(mfile)
  # diff timeseries & metrics! after reading
  s = timeseries[!timeseries %in% metrics$time]
}

# initial state map (= mean pm2.5), 
static$pm25_state = mean(lcs$pm25_smoothed, na.rm=T)

for (i in 1:length(timeseries)){
  
  # read lcs data for current time step
  t = timeseries[i]
  message(t)
  df_now_sf = get_data_now(t, lcs) |> 
    filter(pm25 > 0) |> 
    #filter(pm25 > 0.00001, pm25_lag1 > 0.00001) |> 
    #mutate(across(starts_with("pm25"), log)) |>             ##################
    select(-cell_id, -Land_Cover)
  
  if (nrow(df_now_sf) < 5){
    message("not enough observations. abort!")
    next
  }
  
  # validation data
  val = eea |> filter(time == as.numeric(t))
  
  # add dynamic ERA5/CAMS variables and pm2.5 state
  dynamic_now = dynamic[,,,i, drop=T] |> st_as_stars()
  dynamic_now$pm25_state = static$pm25_state
  dynamic_now_ext = st_extract(dynamic_now, df_now_sf) |> st_drop_geometry()
  df_now_sf = bind_cols(df_now_sf, dynamic_now_ext) |> tidyr::drop_na()
  
  # dynamic PM2.5 interpolation  
  vars = c("pm25_smoothed", "pm25_lag1", "pm25_lag24", "pm25_roll_8h")
  dynamic_pm25 = do.call(
    c, mclapply(vars, nn_interpolation, sf = df_now_sf, 
                  template=static, mc.cores = length(vars)))
  
  # combine static and dynamic covariates
  cov = c(dynamic_pm25, static)

  cov$cams_pm2p5 = dynamic[,,,i, drop=T]$cams_pm2p5
  cov$rh = dynamic[,,,i, drop=T]$rh
  cov$t2m = dynamic[,,,i, drop=T]$t2m
  cov$wd = dynamic[,,,i, drop=T]$wd
  cov$ws = dynamic[,,,i, drop=T]$ws

  # extract covs for validation later
  val_ext = st_extract(cov, val)
  
  # RF model
  fit_classic = randomForest::randomForest(
    x = model.matrix(~ . - 1, data = df_now_sf |> 
                       select(-time, -pm25) |> 
                       st_drop_geometry()),
    y = df_now_sf |> pull("pm25"), #|> log(),            #####################
    importance = TRUE, ntree=100, mtry=3)
  
  fit_classic
  var_imp(fit_classic)
  
  # predict + get residuals
  cov = cov[fa[1,]]
  cov$rf_pred = predict(
    fit_classic, newdata = model.matrix(~ . - 1, data = as.data.frame(cov) |> select(-x,-y))) |> 
    as.vector() #|> exp()                                 ###################
  
  
  #df_now_sf$res = df_now_sf |> pull("pm25") - exp(fit_classic$predicted)
  df_now_sf$res = df_now_sf |> pull("pm25") - fit_classic$predicted
  
  # residual kriging
  cov$res = ord_krig(x=df_now_sf |> filter(!is.na(res)), 
                     var="res", template = cov, 
                     nmin = 5, nmax = 20)
  
  # UK result
  cov$pred = cov$rf_pred + (cov$res) # |> exp())               ################
  #plot(cov[17:19,,] |> merge())
  
  # update state map
  alpha = 0.8  # smoothing factor
  cov["pm25_state"] = alpha * cov["pred"] + (1 - alpha) * cov["pm25_state"]
  
  write_stars(cov["pred"], paste0(map_dir, "/map_", method, "_", as.integer(t), ".tif"))
  write_stars(cov["pm25_state"], paste0(map_dir, "/map_", method, "_state_map_", as.integer(t), ".tif"))
  
  #plot_pm25_pred(cov)
  
  # validation metrics
  # NA if EEA station has no data for this hour
  # predict RF + Kriging only at EEA station location (-> faster!)
  
  
  #val$pred = st_extract(cov["pm25_state"], val) |> pull("pm25_state")
  if (nrow(val) > 0){
    val$rf_pred = predict(
      fit_classic, newdata = model.matrix(~ . - 1, data = val_ext |> st_drop_geometry())) |> 
      as.vector()
    
    val$res = ord_krig(x=df_now_sf |> filter(!is.na(res)), 
                       var="res", template = val_ext, nmin = 5, nmax = 20)
    val$pred = val$rf_pred + val$res
    val$pm25_state = alpha * val$pred + (1 - alpha) * val_ext$pm25_state
  }
  
  metrics = bind_rows(metrics, 
                      list(time = as.character(t), cv = "val", 
                           model_measures(val$raw_pm2_5, val$pred)) |> purrr::flatten(),
                      list(time = as.character(t), cv = "val_state", 
                           model_measures(val$raw_pm2_5, val$pm25_state)) |> purrr::flatten(),
                      list(time = as.character(t), cv = "model", 
                           model_measures(df_now_sf$pm25, fit_classic$predicted)) |> purrr::flatten()
  )
  
  print(metrics[(nrow(metrics)-2):nrow(metrics),])
  saveRDS(metrics, mfile)
  message("------------------------------------")
  gc()
}

list.files(file.path("AQ_maps", method)) |> length()


library(ggplot2)
log_scale = scale_y_continuous(transform = "log", 
                               breaks = c(0,0.01,0.02,0.05,.1,.2,.5,1,2,5,10,20,50,100,150,200)) 


tidyr::pivot_longer(metrics |> select(-MAE, -MBE), 
                    RMSE:R2, names_to = "metric", values_to = "value") |> 
ggplot(aes(x=as.POSIXct(strptime(time, "%Y-%m-%d %H:%M:%S")), 
           y=value, color = cv, pch = metric)) +
  geom_point() +
  labs(title = "PM2.5 - Dynamic State Random Forest - RMSE and R²", x = "Time (hourly)") +
  log_scale


g1 = ggplot(metrics, 
       aes(x=as.POSIXct(strptime(time, "%Y-%m-%d %H:%M:%S")), 
           y=RMSE, color = cv)) +  
  geom_point(alpha = 0.6) +
  labs(title = "PM2.5 - Dynamic State Random Forest - RMSE", x = "Time (hourly)") +
  scale_color_manual("Cross Validation", values = c("deepskyblue", "orange2", "tomato1")) +
  log_scale + theme(legend.position = "bottom") +
  theme_minimal()
g1

g2 = ggplot(metrics |> filter(!cv=="val_state"), 
       aes(x=as.POSIXct(strptime(time, "%Y-%m-%d %H:%M:%S")), 
           y=R2, color = cv)) +  
  geom_point(show.legend = F) +
  labs(title = "PM2.5 - Dynamic State Random Forest - R2", x = "Time (hourly)") +
  scale_color_manual("Cross Validation", values = c("deepskyblue", "orange2", "tomato1")) +
  theme_minimal()

library(patchwork)
g1 + g2 + plot_layout(ncol=1, guides = "collect")
ggsave(file.path("AQ_maps", method, "validation_metrics.png"), width = 5, height = 4)

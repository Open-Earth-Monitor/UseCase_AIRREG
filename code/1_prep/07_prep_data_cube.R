#' ==============================================================================
#' 
#' PREPARE ANALYSIS-READY HOURLY STATION AND COVARIATE DATA
#' 
#' Contents:
#' 1. station loacations
#' 2. station data
#' 3. ERA5 (t2m, rh, ws, wd)
#' 4. CAMS PM2.5
#' 
# ==============================================================================

library(arrow)
library(dplyr)
library(sf)
library(stars)
library(purrr)
library(terra)
library(tictoc)

source("code/1_prep/FILTER/FILTER_QA_QC_Sensor_functions.R")
source("code/airreg_functions.R")

# temporal extent
t1 = "2022-12-31"
t2 = "2024-01-01"

# static covariates
static = readRDS("supplementary/static_predictors.rds") 

# 1. station locations =========================================================

fa = st_read("AQ_data/focus_area.gpkg") # focus area boundary 
stations_fa_sf = st_read("AQ_data/stations_focus_area_30km_buffer.gpkg", quiet=T)[,-2] 

# add lon/lat columns
stations_fa = bind_cols(stations_fa_sf, 
                        st_coordinates(stations_fa_sf) |> 
                          as.data.frame() |> 
                          setNames(c("longitude","latitude"))) |>
  st_transform(st_crs(static))

# attach static covariates
# include cell if -> to avoid kriging artifacts
static["cell_id"] = 1:(dim(static)[1]*dim(static)[2])
ext = st_extract(static, stations_fa) |> st_drop_geometry()
stations_fa = bind_cols(stations_fa, ext) 

st_write(stations_fa, "AQ_data/stations_focus_area_30km_static_covs.gpkg", delete_dsn = T)

# 2. sensor data ===============================================================

from = as.Date(t1) |> as.POSIXct() |> as.integer()
to = as.Date(t2)|> as.POSIXct() |> as.integer()

# OpenSenseMap
osm = open_dataset("AQ_data/hourly/osm_filtered/") |> 
  filter(hour_timestamp >= from, hour_timestamp < to,
         spread > 2, quality_flag == 1) |> 
  inner_join(stations_fa |> st_drop_geometry(), by = c("sensor_id" = "location")) |> 
  rename(time = hour_timestamp, pm25 = raw_pm2_5) |> 
  select(-c(spread, raw_qc, quality_flag)) |> 
  collect() 

# filtered Sensor.Community/PurpleAir
nilu = open_dataset("AQ_data/hourly/nilu/") |> 
  filter(hour_timestamp >= from, hour_timestamp < to, 
         spread > 2) |> 
  inner_join(stations_fa |> st_drop_geometry()) |> 
  mutate(sensor_id = paste(location, sensor_id, sep = "_"),
         pm25 = if_else(is.na(corrected_pm2_5), raw_pm2_5, corrected_pm2_5)) |> 
  collect() |> 
  quality_flag_function(data.column = "raw_pm2_5", qc.column = "raw_qc") |> 
  filter(quality_flag < 3) |> 
  rename(time = hour_timestamp) |> 
  select(all_of(names(osm))) 

# EEA reference
from = as.Date(t1) |> as.POSIXct("CET")
to = as.Date(t2)|> as.POSIXct("CET")
eea = open_dataset("../UseCase_AIRCON/AQ_data/agg_geoparquet/airquality.no2.o3.pm10.pm2p5_1.hourly_pnt_20150101_20231231_eu_epsg.3035_v20240718/") |>
  filter(Start >= from, Start <= to) |>
  collect() |>
  inner_join(stations_fa |> st_drop_geometry(), by = c("Air.Quality.Station.EoI.Code"="location")) |>
  mutate(Air.Quality.Station.EoI.Code = as.character(Air.Quality.Station.EoI.Code),
         hour_timestamp = as.numeric(Start), location = "ref", spread = 9) |>
  filter(Validity_PM2.5 == 1,
         Verification_PM2.5 == 1) |>
  select(time = hour_timestamp, sensor_id = Air.Quality.Station.EoI.Code,
         location, raw_pm2_5 = PM2.5, spread, longitude, latitude)


saveRDS(osm, "AQ_data/aq_data_cube_hourly_2023_opensensemap.rds")
saveRDS(nilu, "AQ_data/aq_data_cube_hourly_2023_sensor_community.rds")
saveRDS(eea, "AQ_data/aq_data_cube_hourly_2023_eea.rds")


# prune data for modelling
lcs = bind_rows(readRDS("AQ_data/aq_data_cube_hourly_2023_sensor_community.rds"),
                readRDS("AQ_data/aq_data_cube_hourly_2023_opensensemap.rds")) |> 
  tidyr::separate(sensor_id, c("s_id", "loc"), sep = "_", remove=F)

stations_fa = st_read("AQ_data/stations_focus_area_30km_static_covs.gpkg") |> 
  rename(s_id = location)

# N observations per sensor
stations_n_obs = group_by(lcs, s_id, loc) |> summarise(n_obs = n())

# keep only stations in lcs (time window of 2023; excludes many stations from earlier)
stations_fa_subset = inner_join(stations_fa, stations_n_obs) 

# N stations per unique grid cell
# should not be more than 1 to avoid kriging artifacts!
cell_counts = stations_fa_subset |> st_drop_geometry() |> count(cell_id)
table(cell_counts$n)

# select one station per unique grid cell 
# based on the number of measurements (n_obs) throughout the time period (e.g., 2023)
# attach counts |> sort by n_obs per grid cell (group) |> keep only station with maximum n_obs
stations_fa_no_duplicates = inner_join(stations_fa_subset, cell_counts) |> 
  arrange(cell_id, -n_obs)|> #filter(loc %in% c("A","B")) |> 
  group_by(cell_id) |> 
  filter(row_number()==1) |> 
  ungroup() |> 
  select(-n)

sum(duplicated(stations_fa_no_duplicates$s_id)) # 0

st_write(stations_fa_no_duplicates, "AQ_data/stations_focus_area_unique_2023.gpkg", delete_dsn = T)

# filter lcs data for unique suitable stations
lcs_no_duplicates = filter(lcs, s_id %in% stations_fa_no_duplicates$s_id & 
                             loc %in% stations_fa_no_duplicates$loc)

lcs_ard = group_by(lcs_no_duplicates, sensor_id) |> 
  arrange(time) |> 
  mutate(
    pm25_smoothed = coalesce(
      zoo::rollapply(pm25, width = 6, FUN = function(x) {
        alpha <- 0.8
        Reduce(function(a, b) alpha * b + (1 - alpha) * a, x, right = TRUE)
      }, fill = NA, align = "right"), pm25),
    pm25_lag1 = coalesce(lag(pm25, 1), pm25),
    pm25_lag24 = coalesce(lag(pm25, 24), pm25),
    pm25_roll_8h = coalesce(zoo::rollmeanr(pm25, 8, fill = NA), pm25)
  ) |> 
  ungroup() |> 
  select(time, longitude, latitude, starts_with("pm25"), all_of(names(static)))

saveRDS(lcs_ard, "AQ_data/lcs_analysis_ready_data_pm25_hourly_2023.rds")

mutate(lcs_ard, 
       # as.single() to convert 8-byte doubles to 4-byte floats
       across(c(starts_with("pm25"), Elevation, Population, contains("Roads"), CLC_Nat), 
              as.single)) |> 
  write_parquet("AQ_data/lcs_analysis_ready_data_pm25_hourly_2023.pq")


# 3. ERA5 ======================================================================

# reference grid
static_rast = as(static["Land_Cover"], "SpatRaster")

from = as.Date(t1) |> as.POSIXct()
to = as.Date(t2)|> as.POSIXct()

era = read_stars(list.files("supplementary/era5_download", ".nc", full.names = T)) |> 
  transmute(
    ws = uv_wind_speed(u10, v10, kmh=T), 
    wd = uv_wind_dir(u10, v10),
    rh = relative_humidity(t2m, d2m),
    t2m = (t2m - units::set_units(273.15, "K"))
  ) |> 
  st_set_dimensions(names = c("lon", "lat", "time")) 

warp_cov = function(name, cov, scale, datatype, threads = 1, 
                    dir = "supplementary", input_crs = "epsg:4326"){
  
  print(name)
  
  if (!inherits(cov, "SpatRaster")){
    t_cov = as(cov[name] |> units::drop_units(), "SpatRaster")
  } else {
    t_cov = cov
  }
  crs(t_cov) = crs(input_crs)
  message("project...")
  t_warp = project(t_cov, static_rast, align = T, method = "bilinear", threads = threads,
                   filename = paste0(dir, "/",name, "_hourly_projected.tif"),
                   gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"), scale = scale, datatype = datatype) 
  gc()
  message("resample...")
  t_warp = resample(t_warp, static_rast, method = "bilinear", threads = threads,
                    filename = paste0(dir, "/",name, "_hourly_resampled.tif"),
                    gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"), scale = scale, datatype = datatype)
  gc()
  return(t_warp)
}

terraOptions(memfrac = 0.8)
#t2m= warp_cov("t2m", era, scale = 0.01, datatype = "INT2S", threads = 2)
#ws = warp_cov("ws", era, scale = 0.01, datatype = "INT2U")
#wd = warp_cov("wd", era, scale = 0.1, datatype = "INT2U")
rh = warp_cov("rh", era, scale = 0.01, datatype = "INT2U", threads = 4)



t = seq(as.POSIXct("2023-01-01T00:00:00", tz = "UTC", format= "%Y-%m-%dT%H"), 
        as.POSIXct("2024-01-01T00:00:00", tz = "UTC", format= "%Y-%m-%dT%H"), by="hour") |> as.numeric()
t = t[-length(t)]
name = "t2m"

st_cov = read_stars("supplementary/era5_t2m_hourly_resampled.tif") |>  
  st_set_dimensions(names = c("x", "y", "time")) |>  
  st_set_dimensions(which ="time", values = t) |> 
  #st_as_stars() |> 
  setNames(name)

saveRDS(st_cov, paste0("supplementary/",name, "_hourly_stars.rds"))



# 4. CAMS ======================================================================

# subset pan-European CAMS to study area, save to NetCDF
rastlist = lapply(list.files("supplementary/cams_download", ".nc", full.names = T),
                  rast, win = ext(7.9, 10, 48.2, 49.2))
t_cams = do.call(c, rastlist)
writeCDF(t_cams, "supplementary/cams_pm2p5_2023.nc", overwrite=T)



tic()
message("project...")
warp_cams = project(t_cams, static_rast, align = T, threads = 6, method = "bilinear", 
                    filename="supplementary/cams_pm25_hourly_2023_projected.tif",
                    gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"), scale = 0.01, datatype = "INT2U")
message("resample...")
warp_cams = resample(warp_cams, static_rast, threads = 6, method = "bilinear", 
                     filename = "supplementary/cams_pm25_hourly_2023_resampled.tif",
                     gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2","TILED=YES"), scale = 0.01, datatype = "INT2U")
print(warp_cams)

st_cams = st_as_stars(warp_cams) |>  
  st_set_dimensions(names = c("x", "y", "time")) |>  
  st_set_dimensions(which ="time", values = cams |> st_get_dimension_values("time")) |> 
  st_as_stars() |> 
  setNames("cams_pm25")

saveRDS(st_cams, "supplementary/cams_pm25_hourly_stars.rds")
toc()

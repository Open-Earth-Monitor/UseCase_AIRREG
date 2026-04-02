library(arrow)
library(dplyr)
library(sf)
library(stars)
library(purrr)
library(terra)
library(tictoc)
library(ggplot2)

source("AIRREG/code/FILTER/FILTER_QA_QC_Sensor_functions.R")
source("AIRREG/code/airreg_functions.R")
source("R/functions.R")

t1 = "2023-02-18"
t2 = "2023-02-26"


# static covariates
static = readRDS("AIRREG/supplementary/static_covariates_buffered.rds") |> 
  setNames(c("Land_Cover", "Elevation", "Population", "Highway", "Major_Roads", "Minor_Roads", "CLC_Nat")) |> 
  mutate(Highway = Highway/1000, Major_Roads = Major_Roads/1000, 
         Minor_Roads = Minor_Roads / 1000)



# station locations
# flag stations within 100 m from other stations - to avoid kriging artifacts

fa = st_read("AIRREG/focus_area.gpkg") # focus area boundary 
stations_fa_sf = st_read("AIRREG/AQ_data/stations_focus_area_30km_buffer.gpkg", quiet=T)[,-2] 
stations_fa_sf$neighbors_within_100m = lengths(st_is_within_distance(stations_fa_sf, dist = 100)) - 1

stations_fa = bind_cols(stations_fa_sf, 
                        st_coordinates(stations_fa_sf) |> 
                          as.data.frame() |> 
                          setNames(c("longitude","latitude"))) |>
  st_transform(st_crs(static))

ext = st_extract(static, stations_fa) |> st_drop_geometry()
stations_fa = bind_cols(stations_fa, ext) |> st_drop_geometry()


# sensor data
t1 = "2023-02-18"
t2 = "2023-02-26"

from = as.Date(t1) |> as.POSIXct() |> as.integer()
to = as.Date(t2)|> as.POSIXct() |> as.integer()

osm = open_dataset("AIRREG/AQ_data/hourly/osm_filtered/") |> 
  filter(hour_timestamp >= from, hour_timestamp < to,
         spread > 2, quality_flag == 1) |> 
  inner_join(stations_fa, by = c("sensor_id" = "location")) |> 
  rename(time = hour_timestamp, pm25 = raw_pm2_5) |> 
  select(-c(spread, raw_qc, quality_flag)) |> 
  collect() 


nilu = open_dataset("AIRREG/AQ_data/hourly/nilu/") |> 
  filter(hour_timestamp >= from, hour_timestamp < to, 
         spread > 2) |> 
  inner_join(stations_fa) |> 
  mutate(sensor_id = paste(location, sensor_id, sep = "_"),
         pm25 = if_else(is.na(corrected_pm2_5), raw_pm2_5, corrected_pm2_5)) |> 
  collect() |> 
  quality_flag_function(data.column = "raw_pm2_5", qc.column = "raw_qc") |> 
  filter(quality_flag < 3) |> 
  rename(time = hour_timestamp) |> 
  select(all_of(names(osm))) 


from = as.Date(t1) |> as.POSIXct()
to = as.Date(t2)|> as.POSIXct()
eea = open_dataset("AQ_data/agg_geoparquet/airquality.no2.o3.pm10.pm2p5_1.hourly_pnt_20150101_20231231_eu_epsg.3035_v20240718/") |>
  filter(Start >= from, Start <= to) |>
  collect() |>
  inner_join(stations_fa |> st_drop_geometry(), by = c("Air.Quality.Station.EoI.Code"="location")) |>
  mutate(Air.Quality.Station.EoI.Code = as.character(Air.Quality.Station.EoI.Code),
         hour_timestamp = as.numeric(Start), location = "ref", spread = 9) |>
  filter(Validity_PM2.5 == 1,
         Verification_PM2.5 == 1) |>
  select(time = hour_timestamp, sensor_id = Air.Quality.Station.EoI.Code,
         location, raw_pm2_5 = PM2.5, spread, longitude, latitude)


# dynamic covariates

from = as.Date(t1) |> as.POSIXct()
to = as.Date(t2)|> as.POSIXct()
era = read_ncdf("AIRREG/spat_stat_cube/era5_cube_raw.nc") |> 
  transmute(ws = as.integer(uv_wind_speed(u10, v10) * 100), 
         wd = as.integer(uv_wind_dir(u10, v10)),
         rh = as.integer(relative_humidity(t2m, d2m)),
         t2m = as.integer((t2m - units::set_units(273.15, "K"))*100)
         ) |> 
  st_set_dimensions(names = c("lon", "lat", "time")) |> 
  filter(between(time, from , to))

cams = read_ncdf("AIRREG/spat_stat_cube/cams.eaq.ira.ENSa.pm2p5.l0.2023-02.nc", 
          var = "pm2p5") |> 
  filter(between(lon, 7, 10),
         between(lat, 47.5, 49.5),
         between(time, from , to)) |> st_as_stars(proxy=F)
cams

identical(cams |> st_get_dimension_values("time"),
          era |> st_get_dimension_values("time"),
          unique(eea$time) |> sort()|> as.POSIXct(),
          unique(nilu$time) |> sort()|> as.POSIXct()) # TRUE



# downsample ERA5
y = as(static, "SpatRaster")

names(era)

write_era = function(n, era){
  print(n)
  tic()
  # era = read_ncdf("AIRREG/spat_stat_cube/era5_cube_raw.nc") |> 
  #   transmute(ws = uv_wind_speed(u10, v10), 
  #             wd = uv_wind_dir(u10, v10),
  #             rh = relative_humidity(t2m, d2m),
  #             t2m = t2m) |> 
  #   st_set_dimensions(names = c("lon", "lat", "time")) |> 
  #   filter(between(time, from , to))
  # 
  t_era = as(era[n] |> units::drop_units(), "SpatRaster")
  message("project...")
  warp_era = project(t_era, y, align = T, method = "bilinear", 
                     filename = paste0("AIRREG/spat_stat_cube/era5_",n, "_hourly_proj.tif")) 
  message("resample...")
  warp_era = resample(warp_era, y, method = "bilinear",
                      filename = paste0("AIRREG/spat_stat_cube/era5_",n, "_hourly_res.tif"))
  
  st_era = st_as_stars(warp_era) |>  
    st_set_dimensions(names = c("x", "y", "time")) |>  
    st_set_dimensions(which ="time", values = era |> st_get_dimension_values("time")) |> 
    st_as_stars() |> 
    setNames(n)

  print(st_era)
  message("write...")
  saveRDS(st_era, paste0("AIRREG/spat_stat_cube/era5_",n, "_hourly_stars.rds"))
  toc()
  message(n, " done.")
  return(st_era)
}

ws = write_era("ws", era)
wd=write_era("wd", era)
rh=write_era("rh", era)
t2m=write_era("t2m", era)

# downsample CAMS
t_cams = as(cams |> units::drop_units(), "SpatRaster")

tic()
message("project...")
warp_cams = project(t_cams, y, align = T, method = "bilinear", filename="AIRREG/spat_stat_cube/cams_pm25_hourly_proj.tif")
message("resample...")
warp_cams = resample(warp_cams, y, method = "bilinear", filename = "AIRREG/spat_stat_cube/cams_pm25_hourly_res.tif")
print(warp_cams)

st_cams = st_as_stars(warp_cams) |>  
  st_set_dimensions(names = c("x", "y", "time")) |>  
  st_set_dimensions(which ="time", values = cams |> st_get_dimension_values("time")) |> 
  st_as_stars() |> 
  setNames("cams_pm25")

saveRDS(st_cams, "AIRREG/spat_stat_cube/cams_pm25_hourly_stars.rds")
toc()



# saveRDS(static, "AIRREG/spat_stat_cube/static_predictors.rds")
# saveRDS(osm, "AIRREG/spat_stat_cube/aq_data_hourly_opensensemap.rds")
# saveRDS(nilu, "AIRREG/spat_stat_cube/aq_data_hourly_sensor_community.rds")
# saveRDS(eea, "AIRREG/spat_stat_cube/aq_data_hourly_eea.rds")

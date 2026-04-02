library(arrow)
library(dplyr)
library(sf)
library(purrr)
library(lubridate)
library(furrr)
plan(multisession(workers = 12))

source("code/1_prep/FILTER/FILTER_QA_QC_Sensor_functions.R")
source("code/1_prep/FILTER/FILTER_QA_QC_Sensor_additional_functions.R")


# DATA =========================================================================

stations_fa = st_read("AQ_data/stations_focus_area_30km_buffer.gpkg") |> 
  select(location)
stations_fa = bind_cols(stations_fa, 
                        st_coordinates(stations_fa) |> 
                          as.data.frame() |> 
                          setNames(c("longitude","latitude"))) |> 
  st_drop_geometry()

#write.csv(stations_fa, "AIRREG/AQ_data/stations_focus_area_30km_buffer.csv")

nilu_files = tibble(
  path = list.files("AQ_data/hourly/nilu", full.names = T),
  date = as.Date(sub(".parquet","", basename(path))))

osm_files = tibble(
  path = list.files("AQ_data/hourly/osm", full.names = T),
  date = as.Date(sub(".parquet","", basename(path))))

dc = "raw_pm2_5"
dt = "hour_timestamp"
time_chunks = seq(as.Date("2018-01-01"), as.Date("2025-01-01")-1, by = "6 month")


for (t in 5:(length(time_chunks))){
  # for ( t in c(13,14)){
  from = time_chunks[t]-16
  to = time_chunks[t+1]+16
  
  message("~~~~~~~~~ Chunk ", t, " ~~~~~~~~~")
  
  osm.files = filter(osm_files, between(date, from, to)) |>
    pull("path") 
  
  if (length(osm.files) > 0){
  osm = open_dataset(osm.files) |> inner_join(stations_fa) |> 
    collect() |> 
    mutate(hour_timestamp = as.integer(hour_timestamp),
           sensor_id = location)
  } else {
    next
  }
  # open only sensor data from the focus area
  # consolidate sensor and location ID
  # derive and filter by quality flag ("high", and "good" for oSM prep)
  nilu = filter(nilu_files, between(date, from, to)) |> 
    pull("path") |> open_dataset() |> 
    inner_join(stations_fa) |> select(all_of(names(osm)), raw_qc) |> 
    mutate(sensor_id = paste(location, sensor_id, sep = "_"))|> collect() |> 
    quality_flag_function(data.column = "raw_pm2_5", qc.column = "raw_qc") |> 
    filter(quality_flag < 3) |> 
    select(-raw_qc, -quality_flag)
  
  from_ms = as.Date(from) |> as.POSIXct()
  to_ms = as.Date(to)|> as.POSIXct()
  
  eea = open_dataset("../UseCase_AIRCON/AQ_data/agg_geoparquet/airquality.no2.o3.pm10.pm2p5_1.hourly_pnt_20150101_20231231_eu_epsg.3035_v20240718/") |> 
    filter(Start >= from_ms, Start <= to_ms) |> 
    collect() |> 
    inner_join(stations_fa |> st_drop_geometry(), by = c("Air.Quality.Station.EoI.Code"="location")) |> 
    mutate(Air.Quality.Station.EoI.Code = as.character(Air.Quality.Station.EoI.Code),
           hour_timestamp = as.numeric(Start), location = "ref", spread = 9) |> 
    filter(Validity_PM2.5 == 1, 
           Verification_PM2.5 == 1) |> 
    select(hour_timestamp, sensor_id = Air.Quality.Station.EoI.Code,
           location, raw_pm2_5 = PM2.5, spread, longitude, latitude)
  
  # SIMPLE QUALITY FLAGS =========================================================
  
  message("Simple Flags oSM")
  q = group_split(osm, sensor_id) |> 
    map(qc_flags, .progress = T) |> 
    list_rbind()  
  
  #table(q$flag_range, useNA = "ifany")
  #table(q$flag_constant_value, useNA = "ifany")
  #table(q$flag_outlier, useNA = "ifany")
  
  message("Simple Flags EEA")
  eea = group_split(eea, sensor_id) |> 
    map(qc_flags, .progress = T) |> 
    list_rbind() |> 
    filter(!if_any(starts_with("flag"), is.na))
  
  
  # SPATIAL QUALITY FLAGS ========================================================
  
  # outlier (neighbors)
  message("Spatial Outliers")
  ids = unique(q$sensor_id)
  q = future_map(ids, .f = function(id){ 
    spatial_outlier_qctest(data = bind_rows(q, nilu, eea), 
                           date.column = dt, data.column = dc, 
                           out.flag.column = "flag_outlier", id.column = "sensor_id", id.target = id, 
                           x.column = "longitude", y.column = "latitude", radii = c(3,30))
  }, .progress = T) |> 
    purrr::list_rbind()
  #table(q$flag_outlier, useNA = "ifany")
  
  
  # spatial correlation
  message("\nSpatial Correlation")
  q = future_map(ids, 
                 .f = function(id){ 
                   spatial_correlation_qctest(data = q, date.column = dt, data.column = dc, #new = F,
                                              id.column = "sensor_id", id.target = id, ndays = 30,
                                              rdata = bind_rows(nilu, eea), rdata.column = dc, rid.column = "sensor_id",
                                              x.column = "longitude", y.column = "latitude")
                 }, .progress = T) |> 
    purrr::list_rbind()
  #table(q$flag_correlation, useNA = "ifany")
  
  
  # spatial similarity
  message("\nSpatial Similarity")
  q = future_map(ids, .f = function(id){ 
    spatial_similarity_qctest(sdata = q, rdata = bind_rows(nilu, eea), id.target = id,
                              date.column = dt, sdata.column = dc, rdata.column = dc,
                              sid.column = "sensor_id", rid.column = "sensor_id",
                              x.column = "longitude", y.column = "latitude",
                              radii = 30, similarity.duration = 168, 
                              compl.threshold = 90, direction='center', param = 'PM2.5')
  }, .progress = T) |> 
    purrr::list_rbind() 
  
  #table(q$flag_similarity, useNA = "ifany")
  
  
  # quality summary
  message("Quality Summary")
  q = quality_flag_function(q, dc) |> 
    filter(hour_timestamp >= as.numeric(as.POSIXct(time_chunks[t], "UTC")), 
           hour_timestamp < as.numeric(as.POSIXct(time_chunks[t+1], "UTC"))) |> 
    mutate(across(starts_with("flag"), ~ tidyr::replace_na(.x, replace = 0)),
           raw_qc = paste0(flag_range, flag_constant_value, flag_outlier, 
                           flag_correlation, flag_similarity),
           date = as.Date(as.POSIXct(hour_timestamp, "UTC"))) |> 
    select(hour_timestamp, sensor_id, raw_pm2_5, spread, raw_qc, quality_flag, date)
  
  message("Writing ", nrow(q), " rows.")
  
  write_dataset(q, path = "AQ_data/hourly/osm_filtered",
                partitioning = "date", hive_style = F)
}

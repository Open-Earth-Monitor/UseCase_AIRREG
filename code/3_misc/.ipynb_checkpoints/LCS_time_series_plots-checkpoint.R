library(lubridate)

t1 = "2022-12-25"
t2 = "2023-12-28"

from = as.Date(t1) |> as.POSIXct() |> as.integer()
to = as.Date(t2)|> as.POSIXct() |> as.integer()


osm = open_dataset("AIRREG/AQ_data/hourly/osm_filtered/") |> 
  filter(hour_timestamp >= from, hour_timestamp < to,
         spread > 2, quality_flag == 1) |> 
  inner_join(stations_fa, by = c("sensor_id" = "location")) |> 
  rename(pm25 = raw_pm2_5) |> 
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
  filter(quality_flag < 2) 

nilu_daily =  mutate(nilu, dt = as.Date(as.POSIXct(hour_timestamp))) |> 
  group_by(sensor_id, dt, Population, Land_Cover, longitude, latitude) |> 
  summarise(pm25 = mean(pm25, na.rm=T))

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
  select(hour_timestamp, sensor_id = Air.Quality.Station.EoI.Code,
         location, pm25 = PM2.5, spread, longitude, latitude)

##################

library(ggplot2)


ggplot(nilu |> 
         # filter(between(hour_timestamp, 
         #                     as.integer(as.POSIXct("2023-02-20")),
         #                     as.integer(as.POSIXct("2023-02-26")))) |> 
         group_by(sensor_id), aes(as.POSIXct(hour_timestamp), pm25)) +
  #geom_line(aes(color=Population), 
            #show.legend = F,
  #          alpha=0.2)
  geom_smooth(span=0.2)

ggplot(nilu_daily |> 
         # filter(between(hour_timestamp, 
         #                     as.integer(as.POSIXct("2023-02-20")),
         #                     as.integer(as.POSIXct("2023-02-26")))) |> 
         group_by(sensor_id), aes(dt, pm25)) +
  geom_boxplot(aes(group=dt),outlier.size = 0.2) +
  scale_y_continuous(limits = c(-5, 100))
            #show.legend = F,
            #alpha=0.2)

library(sf)
library(dplyr)
setwd("~/R/UseCase_AIRCON")

nilu = read.csv("/palma/scratch/tmp/jheisig/aq_reg/data/NILU/Sensor_Location.csv") |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
  mutate(type = as.factor(substr(location, 1, 2)))
  
osm = read.csv("AIRREG/relevant_sensors.csv") |> 
  st_as_sf(coords=c("Lon", "Lat"), crs=4326) |> 
  select(location = Box_ID)|> 
  mutate(type = as.factor("oSM"))

eea = geoarrow::read_geoparquet_sf("AQ_stations/EEA_stations_meta_sf.parquet") |> 
  select(location = Air.Quality.Station.EoI.Code) |> 
  mutate(type = "EEA") |> 
  st_transform(st_crs(osm))
  

all_sensors = rbind(nilu, osm, eea)

plot(all_sensors["type"])
summary(all_sensors$type)

mapview::mapview(all_sensors, zcol = "type")


fa = st_read("AIRREG/focus_area.gpkg")[1,] |> 
  st_buffer(30000) |> 
  st_transform(st_crs(osm))

all_sensors_fa = st_filter(all_sensors, fa[1,]) 
table(all_sensors_fa$type)
write_sf(all_sensors_fa, "AIRREG/AQ_data/stations_focus_area_30km_buffer.gpkg")

plot(all_sensors_fa["type"])
mapview::mapview(all_sensors_fa, zcol = "type", 
                 col.regions = list( "slateblue", "orange", "turquoise4", "black"),
                 cex = 3, alpha=0.4, alpha.regions=0.9) +
  mapview::mapview(fa, alpha.region = 0, color="red", lwd=3)

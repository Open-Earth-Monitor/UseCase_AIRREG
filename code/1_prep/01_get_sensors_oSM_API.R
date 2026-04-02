#===============================================================================
#
# Retrieve PM2.5 Sensors from OpenSenseMap
#
# - use the oSM API to query all boxes
# - find the ones with a PM2.5 sensor
# - add metadata incl. location
#===============================================================================

library(httr)
library(jsonlite)
library(dplyr)


# GET all senseBoxes and their basic meta data through the API
response = GET("https://api.opensensemap.org/boxes")  
status_code(response)
data = content(response, as = "text")  |> fromJSON()


# specify all potential variable names for PM2.5
all_phen = readLines("AIRREG/all_phenomena.txt")
sub_phen = grep("pm2.5", all_phen, ignore.case = T, value = T)

grab_sensor = function(x){
  b = data[x, c("_id", "name", "exposure","createdAt", "model", "currentLocation" )] |> 
    setNames(c("Box_ID", "Box_Name", "Exposure","Created", "Model", "Location"))
  
  # filter for possible phnomenon names
  s = filter(data[x,]$sensors[[1]], title %in% sub_phen) 
  
  # catch the cases where sensor metadata is missing (only important for row binding later)
  nms = c("title", "_id", "sensorType", "unit")
  if (!all(nms %in% names(s))){
    for (n in nms){
      if (!n %in% names(s)){
        s[n] = character(0)
      }
    }
  }
  
  bind_cols(b, select(s, all_of(nms)) |> 
              setNames(c("Phenomenon", "Sensor_ID", "Sensor_Type","Unit"))|> 
              mutate(Sensor_Type = as.character(Sensor_Type))
  )
}



s = purrr::map_dfr(1:nrow(data), grab_sensor)
s$Lon = purrr::map(s$Location$coordinates, 1) |> unlist()
s$Lat = purrr::map(s$Location$coordinates, 2) |> unlist()
s$Location = NULL
s$Phenomenon = "PM2.5"
s$Created = as.Date(s$Created)
s = arrange(s, Created)


# # check sensors
# table(s$Exposure) |> sort()
# table(s$Model) |> sort()
# table(s$Sensor_Type) |> sort()
# table(s$Unit) |> sort()
# plot(s$Lon, s$Lat)
# summary(s$Created)


# filter and export

# filter the right units (µg/m³)
pattern <- "\\(?\\s*[μµumy]g*\\s*/?\\s*m(?:\\^?3|³|ᶟ)\\s*\\)?"

relevant_sensors = filter(s,
       between(Lon, -25, 45),
       between(Lat, 30, 75),
       Exposure == "outdoor",
       str_detect(Model, "luftdaten", negate = T),          # keep only stations that do not appear on Sensor.Community
       str_detect(Unit, pattern)) |> 
  mutate(Box_Name = gsub("[öäü®]", "__", Box_Name),
         Unit = "µg/m³") |> 
  as_tibble()

relevant_sensors$Sensor_Type[str_detect(relevant_sensors$Sensor_Type, "[sS][dDM][sS]?S?.*11")]= "SDS 011"
relevant_sensors$Sensor_Type[str_detect(relevant_sensors$Sensor_Type, "sen5|SEN5")]= "SEN5X"
relevant_sensors$Sensor_Type[str_detect(relevant_sensors$Sensor_Type, "SPS30")] = "SPS30"

write.csv(relevant_sensors, "AIRREG/relevant_sensors.csv")

rm(s, response, data, all_phen, sub_phen, pattern)

# test archive download
# start = min(relevant_sensors$Created)
# rs = relevant_sensors[2,]
# dt = rs$Created
# url = paste0("https://archive.opensensemap.org/",
#              rs$Created, "/", 
#              rs$Box_ID, "-", rs$Box_Name, "/",
#              rs$Sensor_ID, "-", dt, ".csv")
# url
# 
# read.csv(url)



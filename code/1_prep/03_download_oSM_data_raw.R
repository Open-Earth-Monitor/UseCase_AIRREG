#===============================================================================
#
# Download and aggregate measurement data from OpenSenseMap
#
# - read data from URLs
# - write to parquet
#===============================================================================

library(dplyr)
library(readr)
library(purrr)
library(furrr)
library(arrow)

setwd("~/R/UseCase_AIRCON")
link_dir = "AIRREG/AQ_data/links"
pq_dir = "AIRREG/AQ_data/raw"

# URLs from scraping the archive
dates_csvs = list.files(link_dir, pattern = ".csv", full.names = T) |> sort()

# All in one: read one file, aggregate hourly, add Box_ID
read_measurements = function(i, data){
  tryCatch({
    read_csv(data$URL[i],
             col_types = list(col_datetime(), col_double()), 
             show_col_types = FALSE) |> 
      #group_by(Datetime = lubridate::floor_date(createdAt, "hour")) |> 
      #summarize(Value = mean(value, na.rm = TRUE), .groups = "drop") |> 
      mutate(Box_ID = data$Box_ID[i])
  }, error = function(e) return(NULL))
}


# run in parallel for each day
job::job({
  plan(multisession(workers = 10))
       
       for (i in 1:length(dates_csvs)){
         d = sub(".csv","", basename(dates_csvs[[i]]))
         pq_out = file.path(pq_dir, paste0(d, ".parquet"))
         
         if (!file.exists(pq_out)){
           message(d)
           one_day_urls = read_csv(dates_csvs[[i]], show_col_types = F)
           future_map(1:nrow(one_day_urls), read_measurements, data = one_day_urls, .progress = F) |> 
             list_rbind() |> 
             as_arrow_table() |> 
             write_parquet(pq_out)
         }
       }
       plan(sequential)
})

relevant_sensors = read.csv("AIRREG/relevant_sensors.csv")
open_dataset(pq_dir) |> 
  left_join(relevant_sensors |> 
              select(Box_ID, Created, Model, Sensor_Type, Lon, Lat),
            by = "Box_ID") |> 
  collect()

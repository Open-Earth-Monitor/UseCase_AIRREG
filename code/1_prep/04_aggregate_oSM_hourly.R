#===============================================================================
#
# Aggregate measurement data from OpenSenseMap to hourly
#
# - summarize in hourly bins
# - calculate spread (coverage in 20 min bins)
# - write to parquet
#===============================================================================

library(dplyr)
library(furrr)
library(arrow)
library(lubridate)

setwd("~/R/UseCase_AIRCON")
pq_dir = "AIRREG/AQ_data/raw"
out_dir = "AIRREG/AQ_data/hourly/osm"
pq_in = list.files(pq_dir, full.names = T)


# All in one: read one file, aggregate hourly, add Box_ID
aggregate_hourly = function(i, files){
  x = files[i]
  out = file.path(out_dir, basename(x))
  if (!file.exists(out)){
    read_parquet(x, col_types = list(col_datetime(), col_double(), col_character()), 
                 show_col_types = FALSE) |> 
      group_by(location = Box_ID,
               hour_timestamp = floor_date(createdAt, "hour"),
               bin20 = floor_date(createdAt, "20 minutes")) |> 
      # count number of obs per 20-min bin
      summarize(bin_mean = mean(value, na.rm = T),
                n = n(), .groups = "drop") |> 
      # count number of 20-min bins with ≥1 obs
      summarize(raw_pm2_5 = weighted.mean(bin_mean, n, na.rm = T),
                spread = sum(n >= 1), 
                sensor_id = "",
                .by = c(location, hour_timestamp)) |> 
      select(hour_timestamp, sensor_id, location, raw_pm2_5, spread) |> 
      write_parquet(out)
  }
}


# run in parallel for each day
job::job({
  plan(multisession(workers = 10))
  future_map(length(pq_in):1, 
             aggregate_hourly, 
             files = pq_in, 
             .progress = T)
  plan(sequential)
})

pq_done = list.files(pq_dir, full.names = T)



library(arrow)
library(dplyr)
library(lubridate)
library(sf)
library(stringr)
library(purrr)

# stations in the focus area
fa_stations = st_read("AIRREG/AQ_data/stations_focus_area_30km_buffer.gpkg")
table(fa_stations$type)

# SC & PA station CSVs
nilu_csvs = list.files("/palma/scratch/tmp/jheisig/aq_reg/data/NILU/extracted/DEU", full.names = T)
nilu_locations = str_split(basename(nilu_csvs), "_") |> purrr::map(1) |> purrr::list_c()

# ... inside the focus area
nilu_fa_csvs = nilu_csvs[nilu_locations %in% fa_stations$location]

# days of interest
nilu_days = seq(as.Date("2018-01-01"), as.Date("2023-12-31"), by  = 1)

# filter single days from all nilu data and write
rewrite_nilu = function(d, data){
  from = as.POSIXct(d) |> as.integer()
  to = (as.POSIXct(d+1)-1) |> as.integer()
  out = file.path("AIRREG/AQ_data/hourly/nilu", paste0(d, ".parquet"))
  if (!file.exists(out)){
  filter(data, 
         between(hour_timestamp, from, to)) |> 
      write_parquet(out)
  }
}


# map in parallel
job::job({
  library(furrr)
  
  # Split IDs into N chunks (where N = number of workers)
  n_workers = 10
  d_chunks = split(sample(nilu_days), rep(1:n_workers, length.out = length(nilu_days)))
  plan(multisession(workers = n_workers))
  
  # Each worker processes one chunk → only one dataset scan per worker
  result = future_map(d_chunks, function(chunk) {
    
    sch = schema(list(hour_timestamp = int64(),
                      sensor_id = string(),
                      location = string(),
                      raw_pm2_5 = double(),
                      spread = int64(),
                      raw_qc = string(),
                      corrected_pm2_5 = double(),
                      corrected_ci_down =   double(),
                      corrected_ci_up =  double(),
                      correction_improvement = int64(),
                      corrected_qc = string()
    ))
    ds = open_csv_dataset(nilu_fa_csvs, na = c("", "NA", "NULL"), col_types = sch)
    walk(chunk, rewrite_nilu, data = ds)
  }, .progress = T)
  plan(sequential)
})



length(list.files("AIRREG/AQ_data/hourly/nilu")) == length(nilu_days)

# future_map(nilu_days, rewrite_nilu, data = nilu_arrow, .progress = T)
# result = map(nilu_days, rewrite_nilu, data = nilu_arrow, .progress = T)


a = read_parquet("AIRREG/AQ_data/hourly/nilu/2023-04-08.parquet") 
a
summary(a)
h = unique(a$hour_timestamp) |> sort()
as.POSIXct(h, tz="UTC")
length(h)
table(a$hour_timestamp)
table(a$location)


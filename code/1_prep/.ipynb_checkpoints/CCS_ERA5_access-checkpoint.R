library(ecmwfr)

options(keyring_backend="file")
source("code/1_prep/ecmwf_login.R")

wf_set_key(key = key_era, user = user)

dl.dir = "supplementary/era5_download"
if (!dir.exists(dl.dir)) dir.create(dl.dir, recursive = T)

box = c(49.2, 7.9, 48.2, 10.1)

variables = c('10m_u_component_of_wind', 
              '10m_v_component_of_wind',
              '2m_dewpoint_temperature', 
              '2m_temperature')
years = 2023

# create and store requests
for (y in years){
  for(v in variables){
    file = paste0("ERA5_", v, "_hrly_", y,".nc")
    
    request = list(
      date = paste0(y,"-01-01/",y,"-12-31"),
      format = "netcdf",
      variable = v,
      time = paste0(c(paste0(0,0:9), 10:23), ":00"),
      area = box,
      dataset_short_name = "reanalysis-era5-land",
      target = file
    )
    if (!file.exists(file)){
      
        dl = wf_request(
          request  = request,  
          user     = user,
          transfer = F,          # download the file later, not now          
          time_out = 3600*3,
          path     = dl.dir
        )
        print(dl)
        saveRDS(dl, paste0(dl.dir,"/request_",v, "_", y,".rds"))
    }
  }
}

# download in case request has been processed ("completed")
reqs = list.files(dl.dir, pattern="request", full.names = T)
njobs = 1
for (r in reqs) {
  req = readRDS(r)
  req = req$update_status()
  if (req$get_status() == "completed") {
    file = req$get_request()$target
    if (!file.exists(file)){
      message(basename(file), ": start download ", njobs)
      job::job({
        source("code/1_prep/ecmwf_login.R")
        req$.__enclos_env__$private$file = file
        req$.__enclos_env__$private$path = req$.__enclos_env__$private$path |> dirname()
        req$download(verbose = T)
      }, import = "auto", packages = "ecmwfr", title = basename(file))
      njobs = njobs+1
    }
  }
}

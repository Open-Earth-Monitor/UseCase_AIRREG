library(ecmwfr)

options(keyring_backend="file")
source("code/1_prep/ecmwf_login.R")

wf_set_key(key = key_era, user = user)

dl.dir = "supplementary/cams_download"
if (!dir.exists(dl.dir)) dir.create(dl.dir, recursive = T)

v = 'particulate_matter_2.5um'
y = 2023
months = list(1:6, 7:12)

# create and store requests
for (m in months){
  file = paste0("CAMS_", v, "_hrly_", y, "_", min(m), "-", max(m))
  
  request = list(
    variable = v,
    model = "ensemble",
    level = "0",
    type = "validated_reanalysis",
    year = as.character(y),
    month = as.list(as.character(m)),
    dataset_short_name = 'cams-europe-air-quality-reanalyses',
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
    saveRDS(dl, paste0(dl.dir,"/request_",v, "_", y, "_", min(m), "-", max(m), ".rds"))
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
        source("R/ecmwf_login.R")
        req$.__enclos_env__$private$file = file
        req$.__enclos_env__$private$path = req$.__enclos_env__$private$path |> dirname()
        req$download(verbose = T)
      }, import = "auto", packages = "ecmwfr", title = basename(file))
      njobs = njobs+1
    }
  }
}

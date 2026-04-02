library(dplyr)
library(sf)

source("AIRREG/code/FILTER/FILTER_QA_QC_Sensor_functions.R")
source("AIRREG/code/FILTER/FILTER_QA_QC_Sensor_additional_functions.R")
load("AIRREG/code/FILTER/network_sensor_data.rda")
network_sensor_data

nsd = left_join(network_sensor_data$data, network_sensor_data$info, by = "id")



# valid range
nsd_1 = range_validity_qctest(nsd, "PM2.5", upper.bound = 1000, lower.bound = 0)
summary(nsd_1)
table(nsd_1$flag_range, useNA = "ifany")

# constant value
nsd_2 = constant_value_qctest(nsd_1, date.column = "timestamp", data.column = "PM2.5",
                              persist.duration = 8, compl.duration = 1, min.variation = 0.1, 
                              direction = "center", method = "range")
summary(nsd_2$flag_constant_value)
table(nsd_2$flag_constant_value, useNA = "ifany")

# outlier (sensor)
nsd_3 = outlier_qctest(nsd_2, date.column = "timestamp", data.column = "PM2.5",
                       outlier.duration = 360, compl.duration = 1, outlier.threshold = 10,
                       direction='center')
summary(nsd_3)
table(nsd_3$flag_outlier, useNA = "ifany")

# outlier (neighbors)

nsd_4 = purrr::map(unique(nsd_3$id), .f = function(id){ 
  spatial_outlier_qctest(data = nsd_3, date.column = "timestamp", data.column = "PM2.5", 
                         out.flag.column = "flag_outlier", id.column = "id", id.target = id, 
                         x.column = "longitude", y.column = "latitude", radii = c(3,30))
}, .progress = T) |> 
  purrr::list_rbind()
summary(nsd_4)
table(nsd_4$flag_outlier, useNA = "ifany")

# spatial correlation
nsd_5 = purrr::map(unique(nsd_4$id), .f = function(id){ 
  spatial_correlation_qctest(data = nsd_4, date.column = "timestamp", data.column = "PM2.5",
                             id.column = "id", id.target = id, ndays = 30,
                             rdata = NULL, rdata.column = "PM2.5", rid.column = "id",
                             x.column = "longitude", y.column = "latitude")
}, .progress = T) |> 
  purrr::list_rbind()
summary(nsd_5)
table(nsd_5$flag_correlation, useNA = "ifany")



#===============================================================================
#
# Scrape URLs to all PM2.5 measurement files from OpenSenseMap
#
# - crawl through the oSM archive
# - collect URLs to measurement CSVs from predefined sensors
# - run in parallel
#===============================================================================

library(rvest)
library(xml2)
library(stringr)
library(lubridate)
library(dplyr)
library(purrr)
library(furrr)

setwd("~/R/UseCase_AIRCON")
out_dir = "AIRREG/AQ_data/links"

# get relevant SensBoxes (PM2.5 only)
if (file.exists("AIRREG/relevant_sensors.csv")){
  relevant_sensors = read.csv("AIRREG/relevant_sensors.csv")
} else {
  source("AIRREG/code/openSenseMap_API.R")
}

summary(relevant_sensors)

box_ids = relevant_sensors$Box_ID
sensor_ids = relevant_sensors$Sensor_ID

# Function to get links from a given date subdir
get_matching_links <- function(date) {
  
  out_file = file.path(out_dir, paste0(as.Date(date), ".csv"))
  
  if (!file.exists(out_file)){
    base_url <- paste0("https://archive.opensensemap.org/", date, "/")
    
    # Try loading the page
    tryCatch({
      page = read_html(base_url)
      
      # Get all <a href=...> links
      hrefs = page |> html_nodes("a") |> html_attr("href")
      
      # Keep only relevant boxes 
      box_dirs = hrefs[!is.na(hrefs) & !str_detect(hrefs, "^\\./?$")]
      box_dirs = box_dirs[str_detect(box_dirs, paste(box_ids, collapse = "|"))]
      box_dirs = url_absolute(box_dirs, base_url)
      
      # retrieve sensor csv file urls for each box
      csv_list = list()
      for (bd in box_dirs){
        
        # get measurement urls
        bd_hrefs = read_html(bd) |> html_nodes("a") |> html_attr("href")
        
        # pick the right one
        csv_url = bd_hrefs[!is.na(bd_hrefs) & !str_detect(bd_hrefs, "^\\./?$")]
        csv_url = csv_url[str_detect(csv_url, paste(sensor_ids, collapse = "|"))]
        
        # Build full URLs
        csv_url = url_absolute(csv_url, bd)
        
        if (length(csv_url)){
          
          # add meta data
          csv_list = append(csv_list, 
                            list(data.frame("Date" = as.character(as.Date(date)), 
                                            "Box_ID" = str_extract(basename(bd), "^[^-]+"), 
                                            "URL" = csv_url)))
        } 
      }
      
      if (length(csv_list)){
        bind_rows(csv_list) |> readr::write_csv(out_file)
        return(data.frame(Date = as.character(date), Status = "completed"))
        
      } else {
        return(data.frame(Date = as.character(date), Status = "no data of interest"))
      }
      
    }, error = function(e) {
      return(data.frame(Date = as.character(date), Status = "no data"))
    })
  } else {
    return(data.frame(Date = as.character(date), Status = "completed"))
  }
}

# Date range to crawl
start_date = as.Date("2015-01-01")
end_date = as.Date("2024-12-31")

date_range =  seq(start_date, end_date, by = "1 day")

# Get matching sensors and their CSV URL
# run in the background in parallel
job::job({
  plan(multisession, workers=12)
  all_dates = future_map(date_range, get_matching_links, .progress = T) |> list_rbind()
  plan(sequential)
})

table(all_dates$Status)


library(ggplot2)

mutate(all_dates, 
       year = lubridate::year(Date),
       doy = lubridate::yday(Date)) |> 
ggplot(aes(x = doy, y = factor(year), fill = Status)) +
  geom_tile(height = 0.9, alpha=0.5) +
  geom_vline(xintercept = c(lubridate::yday(as.Date(paste0("2023-",2:12,"-01")))),
             linetype = "dashed", color = "slateblue") + 
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Daily Status PM2.5 Data from openSenseMap", fill = "Status",
       y = "Year", x = "Day of Year")





library(ggplot2)
library(dplyr)
setwd("AIRREG/spat_stat_cube")

m1 = readRDS("metrics_RIMM.rds") |> 
  mutate(method = "1. Time-ignorant UK")

m2 = readRDS("metrics_ta_RIMM.rds") |> 
  mutate(method = "2. Time-aware UK")

m3 = readRDS("metrics_RIMM_state_smooth.rds") |> 
  mutate(method = "3. Dynamic State UK")

lcs_avg = readRDS("LowCostSensor_pm25_hourly_Feb23.rds") |> 
  group_by(time) |> 
  summarise(`PM2.5 avg` = mean(pm25, na.rm=T),
            `PM2.5 smooth avg` = mean(pm25_smoothed, na.rm=T)) |> 
  mutate(time = lubridate::as_datetime(time))


metrics = bind_rows(m1,
                    m2,
                    m3) |> 
  mutate(time = lubridate::as_datetime(time),
         cv = if_else(cv == "val", "Validation", "Model")) |> 
  rename(`Cross Validation` = cv)# |>   inner_join(lcs_avg) 

ggplot(metrics, aes(x=time, color = `Cross Validation`)) +
  geom_line(aes(y=RMSE, group = `Cross Validation`)) +
  #geom_line(aes(y=`PM2.5 avg`), color = "red") +
  facet_wrap(~ method, nrow=3) +
  scale_color_viridis_d("Cross Validation",  begin = 0.2, end=0.6) +
  scale_y_continuous(limits = c(0, 35))

ggplot(metrics, aes(x=time, color = `Cross Validation`)) +
  geom_line(aes(y=R2, group = `Cross Validation`)) +
  facet_wrap(~ method, nrow=3) +
  scale_color_viridis_d("Cross Validation", begin = 0.2, end=0.6)


ggplot(metrics, aes(x=time, color = method)) +
  geom_line(aes(y=RMSE)) +
  #geom_line(aes(y=`PM2.5 avg`), color = "red") +
  facet_wrap(~  `Cross Validation`, scales = "free", nrow=2) +
  scale_color_viridis_d("Approach", option = "B", begin = 0.1, end=0.9) +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  labs(x = "")


ggsave("plot_metrics_RMSE.png", width = 9)


ggplot(metrics, aes(x=time, color = method)) +
  geom_line(aes(y=R2)) +
  #geom_line(aes(y=`PM2.5 avg`), color = "red") +
  facet_wrap(~  `Cross Validation`, scales = "free", nrow=2) +
  scale_color_viridis_d("Approach", option = "B", begin = 0.1, end=0.9) +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  labs(x = "")


ggsave("plot_metrics_R2.png", width = 9)


group_by(metrics, `Cross Validation`, method) |> 
  summarise(medRMSE = median(RMSE),
            meanRMSE = mean(RMSE),
            medR2 = median(R2, na.rm=T),
            meanR2 = mean(R2, na.rm=T))


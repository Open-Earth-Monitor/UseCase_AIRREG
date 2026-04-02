library(fpp3)
library(purrr)

lcs2 = bind_rows(readRDS("aq_data_hourly_sensor_community.rds"),
                readRDS("aq_data_hourly_opensensemap.rds")) |> 
  #ilter(neighbors_within_100m == 0) |> 
  group_by(sensor_id) |> 
  arrange(time) |> 
  mutate(
    pm25_smoothed = coalesce(
      zoo::rollapply(pm25, width = 6, FUN = function(x) {
        Reduce(function(a, b) alpha * b + (1 - alpha) * a, x, right = TRUE)
      }, fill = NA, align = "right"), pm25),
    pm25_lag1 = coalesce(lag(pm25, 1), pm25),
    pm25_lag24 = coalesce(lag(pm25, 24), pm25),
    pm25_roll_8h = coalesce(zoo::rollmeanr(pm25, 8, fill = NA), pm25)
  ) |> 
  ungroup() |> 
  select(time, sensor_id, starts_with("pm25"))

lcs2_ts = lcs2 |> st_drop_geometry() |> 
  mutate(time = lubridate::as_datetime(time)) |> 
  tsibble(index = time, key=sensor_id)

fit <- model(lcs2_ts, ETS(pm25_smoothed ~ error("A") + trend("N") + season("N")))
x = map(map(fit$`ETS(pm25_smoothed ~ error("A") + trend("N") + season("N"))`, 1), 1)
x = x[lengths(x)>1]
x = do.call(rbind, x)
x = filter(x, term == "alpha")
hist(x$estimate)
mean(x$estimate)


algeria_economy <- global_economy |>
  filter(Country == "Algeria")
fit <- algeria_economy |>
  model(ETS(Exports ~ error("A") + trend("N") + season("N")))
fc <- fit |> forecast(h = 5)

methods(class=class(fit))
fit$`ETS(Exports ~ error("A") + trend("N") + season("N"))`


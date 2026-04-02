plot_pm25_pred = function(x, pollutant = "PM2.5", layer = "pred", 
                          pts=NULL, ...){
  breaks = c(0, 5, 10, 15, 20, 25, 40)
  
  max_val = dplyr::pull(x, layer) |> max(na.rm = T)
  if (max_val > max(breaks)) breaks[length(breaks)] = max_val+.1
  
  cols = c("darkgreen", "forestgreen", "yellow2", "orange", "red", "darkred")
  titl = paste(pollutant, "Prediction")
  suppressMessages({
    plot(x[layer], col=cols, breaks = breaks, 
         main = titl, key.pos=1, reset = is.null(pts) , ...)
  })
  
  if (!is.null(pts)) plot(pts, add=T, cex=0.4)
}

get_data_now = function(t, dataset){
  df = tidyr::drop_na(dataset) |> 
    dplyr::filter(time == as.integer(t)) |> 
    dplyr::ungroup() |> 
    sf::st_as_sf(coords = c("longitude", "latitude"), crs=4326) |>
    sf::st_transform(3035) |>
    dplyr::group_by(geometry) |>
    dplyr::filter(dplyr::n() == 1) |>
    dplyr::ungroup()
  return(df)
}

relative_humidity <- function(t2m, d2m) {
  # Constants for Magnus formula
  a <- 17.625 
  b <- 243.04 
  const = 273.15 
  
  t2m = units::drop_units(t2m) - const
  d2m = units::drop_units(d2m) - const
  
  # Saturation vapor pressure (es) and actual vapor pressure (e)
  es <- 6.112 * exp((a * t2m) / (b + t2m))
  e  <- 6.112 * exp((a * d2m) / (b + d2m))
  
  # Relative humidity (in %)
  RH <- (e / es) * 100
  
  # Clamp RH between 0 and 100 to handle edge cases
  RH <- pmax(pmin(RH, 100), 0) |> units::set_units("%")

  
  return(RH)
}

uv_wind_speed <- function(u, v, kmh=F) {
  speed <- sqrt(u^2 + v^2)
  
  if (kmh){
    if(units::deparse_unit(speed) %in% c("m s-1", "km h-1")){
      speed = units::set_units(speed, "km/h")
    } else {
      warning("Wind speed unknown, not in [m/s] or [km/h]. Returning values in input unit.")
    }
  }
    
  return(speed)
}

uv_wind_dir <- function(u, v) {
  # atan2 returns angle in radians between -pi and pi; convert to degrees from North
  u = units::drop_units(u)
  v = units::drop_units(v)
  direction <- (atan2(-u, -v) * 180 / pi) %% 360
  direction = units::set_units(direction, "degree")
  return(direction)
}

interpolate_windspeed = function(t, covs){
  
  ws_file = paste0("supplementary/01_daily/wind_speed_daily_", 
                   lubridate::year(as.POSIXct(t)),".nc")
  
  ws = read_stars(ws_file)[,,,lubridate::yday(as.POSIXct(t)), drop=T] |> 
    sf::st_set_crs(sf::st_crs(4326))|> 
    setNames("WindSpeed") |> 
    st_crop(st_bbox(covs) |> st_as_sfc() |> st_transform(4326)) |> 
    st_as_sf(as_points=T) |> 
    st_transform(st_crs(covs))
  
  wsk = ord_krig(x=ws, var="WindSpeed", template = covs, 
                 nmin = 20, nmax = 50) |> 
    st_as_stars() |> setNames("WindSpeed")
  st_dimensions(wsk) = st_dimensions(covs)
  return(c(covs, wsk))
}

interpolate_pm_dynamics = function(targets, points, covs){
  message("Kriging ", paste0(targets, collapse=", "), " in parallel.")
  krig_res = pbmclapply(targets,
                        FUN = function(t){
                          ord_krig(x = points, 
                                   var = t, 
                                   template = covs, 
                                   nmin = 20, nmax = 50)
                        }, mc.cores = length(targets))

  for (i in 1:length(targets)){
    nms = names(covs)
    covs$new = krig_res[[i]] |> as.vector()
    names(covs) = c(nms, targets[i])
  }
  return(covs)
}

model_measures = function(pred, test){
  rmse = r2 = mae = mbe = NA
  if (!any(is.null(pred), is.null(test))){
    rmse = caret::RMSE(as.vector(pred), as.vector(test), na.rm = T) |> round(3)
    r2 = caret::R2(as.vector(pred), as.vector(test), na.rm = T) |> round(3)
    mae = mean(abs(as.vector(pred) - as.vector(test)), na.rm = TRUE) |> round(3)
    mbe = mean(as.vector(pred) - as.vector(test), na.rm = TRUE) |> round(3)
  }
  return(list(RMSE = rmse, R2 = r2, MAE = mae, MBE = mbe))
}

var_imp = function(x){
  df = randomForest::importance(x) |> 
    as.data.frame() |> 
    dplyr::arrange(- `%IncMSE`)
  df$Var = rownames(df)
  rownames(df) = 1:nrow(df)
  dplyr::select(df, 3,1,2)
}

nn_interpolation = function(var, sf, template){
  template = as(template[1], "SpatRaster") 
  dfv = bind_cols(sf |> st_drop_geometry(), 
                  st_coordinates(sf) |> 
                    as.data.frame() |> 
                    setNames(c("x", "y")))
  
  gs <- gstat::gstat(formula=as.formula(paste(var, "~", 1)), 
                     locations = ~x+y,
                     data=dfv, nmax=10, set=list(idp = 0))
  nn <- terra::interpolate(template, gs, debug.level=0)
  return(stars::st_as_stars(nn[[1]]) |> setNames(var))
}


ord_krig = function(x, var, template, nmin = 5, nmax = 10){
  f = as.formula(paste(var, "~", 1))
  v = automap::autofitVariogram(f, x, model = c("Exp", "Sph"))
  k = gstat::krige(f, x, template, v$var_model, nmin = nmin, nmax = nmax, debug.level=0)
  return(k$var1.pred)
}

aviod_Inf = function(x){
  apply(x, 2, na_if, y=Inf) |> 
    apply(2, na_if, y=-Inf)
}




RFGLS_pred_sp = function (RFGLS_out, coords.0, Xtest, h = 1, verbose = FALSE) 
{
  if (missing(coords.0)) {
    stop("error: coords.0 must be specified\n")
  }
  if (!any(is.data.frame(coords.0), is.matrix(coords.0))) {
    stop("error: coords.0 must be a data.frame or matrix\n")
  }
  if (!ncol(coords.0) == 2) {
    stop("error: coords.0 must have two columns\n")
  }
  coords <- RFGLS_out$coords
  func_pred <- RFGLS_predict(RFGLS_out, Xtest, h = h, verbose = verbose)
  func_pred$predicted = aviod_Inf(func_pred$predicted_matrix) |> 
    rowMeans(na.rm = T) #
  
  func_pred_input <- RFGLS_predict(RFGLS_out, RFGLS_out$X, h = h, verbose = verbose)
  func_pred_input$predicted = aviod_Inf(func_pred_input$predicted_matrix) |> 
    rowMeans(na.rm = T)
  
  rfgls_residual <- RFGLS_out$y - func_pred_input$predicted
  est_theta <- BRISC::BRISC_estimation(coords, x = matrix(1, nrow(coords), 1),
                                       y = rfgls_residual, verbose = verbose)
  corr_pred <- BRISC::BRISC_prediction(est_theta, coords.0, 
                                       X.0 = matrix(1, nrow(coords.0), 1), verbose = verbose)
  prediction <- as.vector(corr_pred$prediction + func_pred$predicted)
  
  return(prediction)
}


rfgls_run = function(time, data, static_covs, 
                     pred_t_minus_1, extent = NULL, 
                     mtry = 3, ntree = 10, n.neighbors=15,
                     importance = T, threads = 1){
  
  # filter by time
  message( as.POSIXct(time, "UTC"))
  df_now_sf = get_data_now(time, data)
  
  if (nrow(df_now_sf) < 10){
    message("Not enough measurements at time ", time)
    return(NULL)
  }
  
  # add temporal covariates
  covs = static_covs
  targets = c("smoothed_pm25", "pm25_roll_8h")
  
  #### prediction from previous time step
  if (!is.null(pred_t_minus_1)){
    covs$pm25_lag1 = pred_t_minus_1$pred
    df_now_sf$pm25_lag1 = st_extract(covs["pm25_lag1"], df_now_sf) |> pull(1) 
    message("Using PM2.5 prediction from last timestamp.")
  }
  #### Kriging
  message("Kriging...")
  krig_res = mclapply(targets, 
                      FUN = function(t) ord_krig(x=df_now_sf, var=t, template = covs),
                      mc.cores = length(targets))
  
  for (i in 1:length(targets)){
    nms = names(covs)
    covs$new = krig_res[[i]] |> as.vector()
    names(covs) = c(nms, targets[i])
  }
  
  if (!is.null(extent)) covs = covs[extent]
  
  # RF fit
  
  vars = c("CLC_Nat", "Elevation", "Population", "Highway", "Major_Roads", 
           "Minor_Roads", "smoothed_pm25", "pm25_lag1", "pm25_roll_8h")
  message("RF fit...")
  
  model_mat = model.matrix(~ . - 1, data = df_now_sf |> 
                             select(all_of(vars)) |> 
                             st_drop_geometry())
  suppressWarnings({
    fit = RFGLS_estimate_spatial(
      coords = st_coordinates(df_now_sf),
      y = df_now_sf$pm25,
      X = model_mat,
      mtry = mtry,
      ntree = ntree,
      n_omp = threads,
      n.neighbors = n.neighbors,
      cov.model = "exponential",
      param_estimate = T,
      verbose = F
    )
    fit$predicted = aviod_Inf(fit$predicted_matrix) |> rowMeans(na.rm = T)
    
    m = data.frame(t = as.POSIXct(time, "UTC"), 
                   n = nrow(df_now_sf),
                   mtry = mtry, ntree=ntree, nn=n.neighbors) |> 
      cbind(as.data.frame(model_measures(fit$predicted, fit$y)))
    print(m)
    
    message("Variable importance...")
    var_imp = NULL
    if (importance){
      var_imp = randomForest::randomForest(
        x = model_mat,
        y = df_now_sf$pm25, 
        importance = TRUE, ntree=ntree, mtry=mtry) |> 
        randomForest::importance()
      var_imp = cbind(t = as.POSIXct(time, "UTC"), var=rownames(var_imp), var_imp) |> 
        as.data.frame(row.names = F)
    }
    
    message("RF prediction...")
    covs$pred = RFGLS_pred_sp(RFGLS_out = fit, 
                              coords.0 = as.matrix(st_coordinates(covs)),
                              Xtest = model.matrix(~ . - 1, 
                                                   data = as.data.frame(covs) |> 
                                                     select(all_of(vars))))
    covs$pred[covs$pred < 0] = 0
  })
  return(list(fit=fit, metrics=m, var_imp=var_imp, pred = covs["pred"]))
}

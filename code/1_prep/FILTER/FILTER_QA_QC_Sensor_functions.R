# range validity test: ensures that observations fall within acceptable range limits
# requires lower and upper limit
# NA: no response -> NA (Not Available) or NaN (Not a Number) 
# 0: not valid -> if the measurement is outside the measurement range 
# 1: valid -> if the measurement is within the measurement range
range_validity_qctest <- function(data, 
                                  data.column, 
                                  upper.bound, 
                                  lower.bound, 
                                  param = 'PM2.5',
                                  new = TRUE) {
  # argument flags
  stopifnot(is.numeric(upper.bound), is.numeric(lower.bound), !is.null(upper.bound), !is.null(lower.bound))
  stopifnot(upper.bound > lower.bound)
  stopifnot(is.character(data.column), data.column %in% colnames(data))
  #
  obs.data <- pull(data,data.column)
  output.data <- data %>% 
    mutate(flag_range = if_else(!!sym(data.column) >= lower.bound  & !!sym(data.column) <= upper.bound, true = 1, false= 0, missing = NA))
  if (new  == TRUE) {
    output.data <- output.data %>% 
      mutate(!!sym(paste0(param,'_range')) := ifelse(flag_range == 0 | is.na(flag_range), NA, obs.data))
  }
  output.data <- output.data %>% as_tibble()
  return(output.data)
}
#
#
# persistence test: detects unchanged behavior over pre+defined time window (similar consecutive values)
# detects censored data, censoring of low measurements (i.e. consecutive zeros over short time periods), 
# detects spurius data caused by sensor's malfunction (non-zero or long-lasting zeros)
# requires time window (persist.duration) and threshold (min.variation)
# two options: min-max range or standard deviation
# NA: no response -> the value is NA (Not Available) or NaN (Not a Number)
# 0: not adequate data -> data completeness within the temporal window is below the completeness threshold (user defined) 
# 1: constant values -> the calculated statistic is lower than the persistence threshold (user defined) and 
#    data completeness within the temporal window is above the completeness threshold (user defined)
# 2: non constant values -> the calculated statistic is higher than the persistence threshold (user defined) and 
#    data completeness within the temporal window is above the completeness threshold (user defined)
constant_value_qctest <- function(data, 
                                  date.column, 
                                  data.column,
                                  persist.duration, 
                                  compl.duration, 
                                  min.variation, 
                                  direction = c('left','center','right'), 
                                  method=c('std','range'),
                                  param = 'PM2.5',
                                  new = TRUE) {
  # argument flags
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(data.column), data.column %in% colnames(data))
    stopifnot(persist.duration > 0, !is.null(persist.duration))
  stopifnot(compl.duration >= 0, !is.null(compl.duration))
  stopifnot(min.variation >= 0)
  stopifnot(direction %in% c('left','center','right'))
  stopifnot(method %in% c('range','std'))
  
  # constant_value check
  compl.duration <- ifelse(compl.duration < 1, round(compl.duration*persist.duration), compl.duration)
  #obs.data <- pull(data,data.column)
  obs.data = select(data, all_of(c(date.column, data.column)))
  ind.nona <- which(!is.na(pull(obs.data,data.column))) 
  
  if (length(ind.nona) == 0) {
    output.data <- data %>% mutate(flag_constant_value = NA)
  } else {
    start.nona <- pull(data,date.column)[ind.nona[1]]
    end.nona <- pull(data,date.column)[ind.nona[length(ind.nona)]]
    data.sub <- data %>% filter(!!sym(date.column) >= start.nona & !!sym(date.column) <= end.nona)
    
    if (method == 'range') {
      statistic <- zoo::rollapply(pull(data.sub,data.column), persist.duration, 
                                  FUN=function(x) {
                                    max.x <- ifelse(sum(is.na(x)) == length(x),NA,max(x,na.rm=TRUE))
                                    min.x <- ifelse(sum(is.na(x)) == length(x),NA,min(x,na.rm=TRUE))
                                    range.x <- ifelse(is.na(max.x) | is.na(min.x) | is.infinite(max.x) | is.infinite(min.x) , NA, round(max.x - min.x,1))
                                    return(range.x)}, 
                                  by.column= FALSE,fill = NA,align=direction)

    } else if (method == 'std') {
      statistic <- zoo::rollapply(pull(data.sub, data.column), persist.duration, 
                                  FUN=function(x) round(sd(x,na.rm=TRUE),1), 
                                  by.column= FALSE, fill = NA, align=direction)
    }
    
    nobs.statistic <- zoo::rollapply(pull(data.sub,data.column), persist.duration, 
                                     FUN=function(x) { 
                                       nobs <- sum(!is.na(x))
                                       return(nobs)}, 
                                     by.column= FALSE,fill = NA,align=direction)
    output.data <- mutate(data.sub, flag_constant_value = 
                            ifelse(nobs.statistic < compl.duration, 0, 
                                   ifelse(statistic <= min.variation & nobs.statistic >= compl.duration, 1,
                                          ifelse(statistic > min.variation & nobs.statistic >= compl.duration, 2, NA))))
    
    output.na.startend <- data %>% 
      filter(!!sym(date.column) < start.nona | !!sym(date.column) > end.nona)
    
    output.data <- bind_rows(output.na.startend, output.data) %>% 
      as_tibble() |> 
      arrange(!!sym(date.column)) |> 
      mutate(flag_constant_value = ifelse(is.na(!!sym(data.column)), NA, flag_constant_value))
  }
  if (new == TRUE) {
    output.data = output.data %>% 
      mutate(!!sym(paste0(param, '_constant_value')) := 
               ifelse(flag_constant_value == 0 | 
                        flag_constant_value == 1 | 
                        is.na(flag_constant_value), NA, !!sym(data.column)))
  }

  return(as_tibble(output.data))
}
#
# Outlier test
# outliers are defined as measurements more or less than 10 times the median absolute deviations (MAD) from the median within 
# a pre-defined temporal window
# NA: no response -> the value is NA (Not Available) or NaN (Not a Number) 
# 0: not adequate data -> data completeness within the temporal window is below the completeness threshold (user defined) 
# 1: outlier value -> the calculated statistic is higher than the outlier threshold (user defined) and data completeness within
# the temporal window is above the completeness threshold (user defined)
# 3: non outlier value -> the calculated statistic is lower than the outlier threshold (user defined) and 
# data completeness within the temporal window is above the completeness threshold (user defined)
outlier_qctest <- function(data, 
                           date.column, 
                           data.column, 
                           outlier.duration, 
                           compl.duration, 
                           outlier.threshold,
                           direction=c('left','center','right'),
                           param = 'PM2.5',
                           new = TRUE) {
  # argument flags
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(data.column), data.column %in% colnames(data))
  stopifnot(outlier.duration > 0, !is.null(outlier.duration))
  stopifnot(compl.duration >= 0, !is.null(compl.duration))
  stopifnot(direction %in% c('left','center','right'))
  stopifnot(outlier.threshold >= 0)
  # outlier check
  compl.duration <- ifelse(compl.duration < 1, round(compl.duration*outlier.duration), compl.duration)
  obs.data <- pull(data,data.column)
  ind.nona <- which(!is.na(obs.data)) 
  if (length(ind.nona) == 0) {
    output.data <- data %>% 
      mutate(flag_outlier = NA)
  } else {
    start.nona <- pull(data,date.column)[ind.nona[1]]
    end.nona <- pull(data,date.column)[ind.nona[length(ind.nona)]]
    data.sub <- data %>% 
      filter(!!sym(date.column) >= start.nona & !!sym(date.column) <= end.nona)
    #
    roll.diff.med <- zoo::rollapply(pull(data.sub,data.column), outlier.duration, FUN=function(x) {
      median.x <- ifelse(sum(is.na(x)) == length(x),NA,round(median(x,na.rm=TRUE),1))
      z.x <- ifelse(is.na(median.x) | is.infinite(median.x) , NA, abs(round(x[length(x)] - median.x,1)))
      return(z.x)}, by.column= FALSE,fill = NA,align=direction)
    #
    roll.mad <- zoo::rollapply(pull(data.sub,data.column), outlier.duration, FUN=function(x) {
      mad.x <- ifelse(sum(is.na(x)) == length(x),NA,round(mad(x,na.rm=TRUE,constant = 1),1))
      mad.x <- ifelse(is.na(mad.x) | is.infinite(mad.x), NA, mad.x)
      return(mad.x)}, by.column= FALSE,fill = NA,align=direction)
    #
    roll.nobs <- zoo::rollapply(pull(data.sub,data.column), outlier.duration, FUN=function(x) { 
      nobs <- sum(!is.na(x))
      return(nobs)}, by.column= FALSE,fill = NA,align=direction)
    output.data <- data.sub %>% 
      mutate(flag_outlier = ifelse(roll.nobs < compl.duration, 0, 
                                          ifelse(roll.diff.med >= outlier.threshold*roll.mad & roll.nobs >= compl.duration, 1,
                                                 ifelse(roll.diff.med < outlier.threshold*roll.mad & roll.nobs >= compl.duration, 3, NA))))   ## adjusted 2 to 3
    #
    output.na.startend <- data %>% 
      filter((!!sym(date.column) < start.nona) | (!!sym(date.column) > end.nona))
    output.data <- bind_rows(output.na.startend,output.data) %>% 
      as_tibble() 
    output.data <- output.data %>% 
      arrange(!!sym(date.column))
    #
    output.data <- output.data %>% 
      mutate(flag_outlier = ifelse(is.na(!!sym(data.column)),NA,flag_outlier))
  }
  if (new == TRUE) {
    output.data <- output.data %>% 
      mutate(!!sym(paste0(param, '_outlier')) := ifelse(flag_outlier == 0 | flag_outlier == 1 | is.na(flag_outlier), NA, obs.data))
  }
  output.data <- output.data %>% 
    as_tibble()
  #
  return(output.data)
}
#
#
# Spatial outlier test
# Outliers at the same timestamps from sensors in the neighborhood spatial scale (~250 m radii centered to sensor locations) 
# are specified events and not outliers 
# 1: Outlier -> the target sensor does not show outliers at the same timestamps with the neighboring sensors 
# 2: Not an outlier -> target and neighboring sensors report outliers at the same timestamp: specific events
spatial_outlier_qctest <- function(data,
                                   date.column,
                                   data.column,
                                   out.flag.column,
                                   id.column,
                                   id.target,
                                   x.column,
                                   y.column,
                                   radii = c(2,25),
                                   param = 'PM2.5',
                                   new = TRUE) {
  # argument flags
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(data.column), data.column %in% colnames(data))
  stopifnot(is.character(out.flag.column), out.flag.column %in% colnames(data))
  stopifnot(is.character(id.column), id.column %in% colnames(data))
  stopifnot(!is.null(id.target))
  stopifnot(is.character(x.column),x.column %in% colnames(data))
  stopifnot(is.character(y.column),y.column %in% colnames(data))
  stopifnot(is.numeric(radii), length(radii) > 1, !is.null(radii))
  #
  # maximum radii in m
  radii <- radii*1000
  
  # target location
  target.data <- filter(data, 
                        !!sym(id.column) %in% id.target,
                        !!sym(out.flag.column) == 1) %>% 
    select(!!sym(date.column),!!sym(data.column),
           !!sym(id.column),!!sym(out.flag.column))   # %>% rename_at(c(data.column,id.column), list( ~paste0(., '_target')))
  
  if (nrow(target.data) == 0) {
    output.data <- data %>% filter(!!sym(id.column) %in% id.target)
  } else {
    target.info <- data %>% filter(!!sym(id.column) %in% id.target) %>% 
      select(!!sym(id.column),!!sym(x.column),!!sym(y.column)) %>% unique()
    target <- sf::st_as_sf(target.info, crs = 4326L, coords = c(x.column,y.column))
    
    # neighboring locations
    # neighbor data 
    neighbors.data <- filter(data, 
                             !(!!sym(id.column) %in% id.target),
                             !!sym(date.column) %in% pull(target.data, date.column))
    if (nrow(neighbors.data) == 0) {
      output.data <- data %>% filter(!!sym(id.column) %in% id.target) 
    } else {
      neighbors.data <- filter(data, 
                               !(!!sym(id.column) %in% id.target),
                               !!sym(date.column) %in% pull(target.data, date.column),
                               !!sym(out.flag.column) == 1) %>%
        select(!!sym(date.column),!!sym(data.column),!!sym(id.column),!!sym(out.flag.column)) 
      
      if (nrow(neighbors.data) == 0) {
        output.data <- data %>% filter(!!sym(id.column) %in% id.target)
      } else  {
        neighbors.info <- data %>% 
          filter(!!sym(id.column) %in% unique(pull(neighbors.data, id.column))) %>%
          select(!!sym(id.column),!!sym(x.column),!!sym(y.column)) %>% unique()
        neighbors <- sf::st_as_sf(neighbors.info, crs = 4326L, coords = c(x.column,y.column))
        
        # calculation of separation distance and selection of neighbors
        distances <- round(sf::st_distance(target,neighbors),1)
        neighbors.info <- neighbors.info %>% 
          mutate(dij = as.vector(distances)) %>% 
          arrange(dij) %>% 
          filter(dij <= radii[2])
        if (nrow(neighbors.info) == 0) {
          output.data <- data %>% filter(!!sym(id.column) %in% id.target)
        } else {
          # calculation of the spatial outlier condition
          # local scale: <= 2 km 
          # not-local scale: > 2 km
          neighbors.data <- neighbors.data %>% 
            filter(!!sym(id.column) %in% pull(neighbors.info, id.column))
          neighbors.data <- plyr::join(neighbors.data, neighbors.info %>% 
                                         select(c(!!sym(paste0(id.column)),dij)), 
                                       by = id.column) %>% as_tibble()
          
          target.neighbors.data <- plyr::join(target.data %>% select(!!sym(date.column), 
                                                                     !!sym(id.column)) %>% rename(!!sym(paste0(id.column,'_target')) := !!sym(id.column)), 
                                              neighbors.data %>% select(!!sym(date.column), !!sym(id.column),'dij') %>% rename(!!sym(paste0(id.column,'_neigh')) := !!sym(id.column)),
                                              by=date.column) %>% 
            as_tibble() %>% tidyr::drop_na() |> 
            mutate(radii_cut = ifelse(dij <= radii[1], 'local', 'not_local'))
          
          spatial.outlier.statistic <- group_by(target.neighbors.data, !!sym(date.column)) %>% 
            summarise(local = length(which(radii_cut == 'local')), 
                             not_local = length(which(radii_cut == 'not_local'))) %>% 
            mutate(!!sym(paste0(out.flag.column,'_spatial')) := ifelse(local >= 1 | not_local >= 2, 2, 1))
          
          output.data <- filter(data, !!sym(id.column) %in% id.target) %>% 
            mutate(!!sym(out.flag.column) := ifelse(!!sym(date.column) %in% pull(spatial.outlier.statistic, date.column), 
                                                           pull(spatial.outlier.statistic, paste0(out.flag.column,'_spatial')), 
                                                           !!sym(out.flag.column)))
        }
      }
    }
  }
  output.data <- output.data %>% 
    mutate(!!sym(out.flag.column) :=  ifelse(is.na(!!sym(data.column)), NA, !!sym(out.flag.column)))
  if (new == TRUE) {
    output.data <- output.data %>% 
      mutate(!!sym(paste0(param, '_outlier')) := ifelse(!!sym(out.flag.column) == 0 | !!sym(out.flag.column) == 1 | is.na(!!sym(out.flag.column)), NA, !!sym(data.column)))
  }
  output.data <- output.data %>% 
    as_tibble()
  #
  return(output.data)
}
#
#
# spatial correlation between sensors
# 0: not adequate data -> data completeness within the temporal window (user defined) is below the completeness threshold (user defined) 
# 1: Not enough neighbors -> No neighbor within 3 km and only one neighbor within 30 km
# 2: No Correlation -> The target location is not correlated with the neighbors because of malfunction issues of lack of correlation due to spatial variability
# 3: Correlation -> The target location is correlated with the nearby locations
# 
spatial_correlation_qctest <- function(data,
                                        rdata,
                                        date.column, 
                                        data.column,
                                        rdata.column,
                                        id.column, 
                                        rid.column,
                                        id.target, 
                                        x.column, 
                                        y.column,
                                        radii = c(3,30),
                                        minneigh = 2, 
                                        ndays = 30,
                                        compl.threshold = 90,
                                        param = 'PM2.5',
                                        new = TRUE) {
  # argument flags
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(data.column), data.column %in% colnames(data))
  stopifnot(is.character(id.column), id.column %in% colnames(data))
  #stopifnot(is.character(id.target))
  stopifnot(is.character(x.column),x.column %in% colnames(data))
  stopifnot(is.character(y.column),y.column %in% colnames(data))
  stopifnot(is.numeric(radii), !is.null(radii))
  stopifnot(is.numeric(minneigh), !is.null(minneigh))
  stopifnot(is.numeric(ndays), ndays > 0, !is.null(ndays))
  stopifnot(is.numeric(compl.threshold), !is.null(compl.threshold))
  #
  # maximum radii in m
  radii <- radii*1000
  compl.threshold <- ifelse(compl.threshold < 1, round(compl.threshold*ndays*24), compl.threshold)
  
  # reference data
  if (!is.null(rdata)) {
    stopifnot(!is.null(rdata.column), is.character(rdata.column), rdata.column %in% colnames(rdata))
    stopifnot(!is.null(date.column), is.character(date.column), date.column %in% colnames(rdata))
    stopifnot(!is.null(rid.column), is.character(rid.column), rid.column %in% colnames(rdata))
    
    rdata <- dplyr::rename(rdata, !!sym(id.column) := !!sym(rid.column),
                           !!sym(data.column) := !!sym(rdata.column)) %>% 
      dplyr::select(all_of(c(date.column, id.column, x.column, y.column, data.column)))
    
    data <- bind_rows(data, rdata) %>% as_tibble()
  }
  
  # mask previously flagged obs
  mutate(data, )
  
  # target location
  target.data <- data %>% filter(!!sym(id.column) %in% id.target) %>% 
    dplyr::select(!!sym(date.column),!!sym(data.column),!!sym(id.column))
  target.info <- data %>%
    filter(!!sym(id.column) %in% id.target) %>%
    dplyr::select(!!sym(id.column),!!sym(x.column),!!sym(y.column)) %>% unique()
  target <- st_as_sf(target.info, crs = 4326L, coords = c(x.column,y.column))
  
  # split the target data into temporal intervals of 30 days (starting at 23:00)
  days = pull(target.data, !!sym(date.column)) |> 
    as.POSIXct(tz="UTC") |> 
    lubridate::as_date() |> 
    unique() 
  days.intervals = paste(days, '23:00:00') |> lubridate::ymd_hms()
  
  target.temporal.windows <- runner::runner(
    x = data.frame(target.data),
    f = function(x) x %>% as_tibble(), 
    k = paste(ndays,'days'),
    idx = date.column,
    at = days.intervals,
    simplify = FALSE
  )
  
  days.spatial.correlation.flag <- rep(NA, length.out=length(target.temporal.windows))
  for (i in 1:length(target.temporal.windows)) {
    target.temporal.tbl <- target.temporal.windows[[i]]
    
    if (sum(!is.na(pull(target.temporal.tbl,data.column))) < compl.threshold) {
      days.spatial.correlation.flag[i] <- 0
    } else {
      # neighbor data 
      # omit neighbor data outside the temporal window
      neighbors.temporal.data <- data %>% 
        filter(!(!!sym(id.column) %in% id.target),
               !!sym(date.column) >= min(pull(target.temporal.tbl, date.column)) & 
                 !!sym(date.column) <= max(pull(target.temporal.tbl, date.column)))
      
      # find neighbors with sufficient paired measurements
      neighbors.temporal.data.list <- neighbors.temporal.data %>% 
        group_by(!!sym(id.column)) %>% 
        group_split() |> 
        lapply(function(x) {
          if (sum(!is.na(pull(x, data.column))) < compl.threshold) {
            return(NULL)
          } else {
            return(x)
          }
        })
      neighbors.temporal.data.list <- rlist::list.clean(neighbors.temporal.data.list, fun = is.null)
      neighbors.temporal.data <- bind_rows(neighbors.temporal.data.list) %>% as_tibble() 
      
      if (nrow(neighbors.temporal.data) == 0) {
        days.spatial.correlation.flag[i] <- 0
      } else {
        
        neighbors.temporal.tbl <- neighbors.temporal.data %>% 
          dplyr::select(!!sym(date.column),!!sym(data.column),!!sym(id.column))
        neighbors.temporal.info <- neighbors.temporal.data %>% 
          dplyr::select(!!sym(id.column),!!sym(x.column),!!sym(y.column)) %>% 
          unique()
        neighbors <- sf::st_as_sf(neighbors.temporal.info, crs = 4326L, coords = c(x.column,y.column))
        
        # calculation of separation distance and selection of neighbors
        distances <- round(st_distance(target,neighbors), 1)
        neighbors.temporal.info <- neighbors.temporal.info %>% 
          dplyr::mutate(dij = as.vector(distances)) %>% 
          dplyr::arrange(dij) %>% 
          dplyr::filter(dij < radii[2])
        
        # neighboring and target data in matrix format
        neighbors.temporal.tbl <- neighbors.temporal.tbl %>% 
          filter(!!sym(id.column) %in% pull(neighbors.temporal.info, id.column))
        
        neighbors.temporal.tbl.mat <- tidyr::pivot_wider(
          neighbors.temporal.tbl %>% 
            dplyr::select(all_of(c(date.column, id.column, data.column))), 
          names_from = !!sym(id.column), 
          values_from = !!sym(data.column)) %>% as_tibble()
        
        neighbors.temporal.tbl.mat <- plyr::join(
          target.temporal.tbl %>% 
            dplyr::select(all_of(c(date.column, data.column))), 
          neighbors.temporal.tbl.mat, by=date.column) %>% as_tibble() %>% 
          dplyr::select(-all_of(c(date.column, data.column)))
        
        target.temporal.tbl.mat <- replicate(ncol(neighbors.temporal.tbl.mat), 
                                             pull(target.temporal.tbl, data.column))
        
        # calculation of pairwise correlation and completeness
        completeness.tbl.mat <- apply(target.temporal.tbl.mat + neighbors.temporal.tbl.mat,
                                      2, function(x) sum(!is.na(x)))
        completeness.stat <- tibble(!!sym(id.column) := names(completeness.tbl.mat), 
                                    completeness := as.vector(completeness.tbl.mat))
        pairwise.correlation <- mapply(cor, 
                                       neighbors.temporal.tbl.mat %>% as.data.frame(), 
                                       target.temporal.tbl.mat %>% as.data.frame(), 
                                       method='pearson', use='pairwise.complete.obs')
        correlation.stat <- tibble(!!sym(id.column) := names(pairwise.correlation), 
                                   correlation := as.vector(pairwise.correlation))
        correlation.tbl <- plyr::join(correlation.stat, completeness.stat, by=id.column) %>% as_tibble()
        
        # calculation of the expected correlation
        expected.correlation.tbl <- neighbors.temporal.tbl %>% 
          dplyr::mutate(Season = quarter(as.POSIXct(!!sym(date.column), "UTC"), fiscal_start = 12)) %>% 
          dplyr::select(all_of(c('Season',id.column))) %>% 
          dplyr::group_by(!!sym(id.column)) %>% 
          dplyr::filter(Season == min(DescTools::Mode(Season))) %>% distinct()
        
        expected.correlation.tbl <- plyr::join(expected.correlation.tbl, 
                                               neighbors.temporal.info %>% 
                                                 dplyr::select(all_of(c(id.column, 'dij'))), 
                                               by=id.column) %>% as_tibble()
        expected.correlation.tbl <- expected.correlation.tbl %>% 
          dplyr::mutate(
            expected_correlation = expected_correlation_fun(season = Season, 
                                                            distance = dij/1000, 
                                                            relax.factor = 0.25)) 
        
        # comparisons
        comparison.tbl <- plyr::join(correlation.tbl, 
                                     expected.correlation.tbl %>% 
                                       dplyr::select(-Season), by = id.column) |> 
          filter(completeness > compl.threshold)
        
        if (nrow(comparison.tbl) == 0) {
          days.spatial.correlation.flag[i] <- 1
        } else {
          comparison.tbl <- comparison.tbl %>% 
            filter(correlation > expected_correlation)
          
          days.spatial.correlation.flag[i] <- ifelse(
            (sum(pull(comparison.tbl, 'dij') <= radii[1], na.rm = TRUE) >= 1) | 
              (sum(pull(comparison.tbl, 'dij') > radii[1], na.rm = TRUE) >= minneigh), 3, 2)
        }
      }
    }
  }
  # daily spatial correlation flags 
  # addition of correlation flag to the sensor data
  daily.correlation.flag <- data.frame(Date = days, 
                                       flag_correlation = days.spatial.correlation.flag)
  output.data <- data %>% filter(!!sym(id.column) %in% id.target) %>% 
    dplyr::mutate(Date = as_date(as.POSIXct(!!sym(date.column), "UTC"))) |> 
    plyr::join(daily.correlation.flag, by='Date') %>% 
    as_tibble() |> dplyr::select(-Date) |> 
    dplyr::mutate(flag_correlation := ifelse(is.na(!!sym(data.column)), NA, flag_correlation))
  
  if (new == TRUE) {
    output.data <- output.data %>% 
      mutate(!!sym(paste0(param,'_correlation')) := 
               ifelse(flag_correlation != 3 | is.na(flag_correlation), 
                      NA, pull(target.data, data.column)))
  }
  
  return(as_tibble(output.data))
}


#
#
# Spatial consistency check
# Similarity of sensor measurements to reference sites
# NA: no response -> the value is NA (Not Available) or NaN (Not a Number)
# 0: No nearby sensor within the neighborhood (~250 m)
# 1: not adequate data -> data completeness within the temporal window (user defined) is below the completeness threshold (user defined) 
# 2: not consistent data -> not consistent measurements over the temporal window (user defined) between the target and the network of sensors within the neighborhood (~250 m radii cenetered to the sensor locations)
# 3: consistent data: -> -> not consistent measurements over the temporal window (user defined) between the target and the network of sensors within the neighborhood (~250 m radii cenetered to the sensor locations)
# (non-) consistency is determined through: Coefficient.of-Divergence (CoD) and Absolute Difference (AD)
#
spatial_similarity_qctest <- function(sdata,
                                      rdata,
                                      date.column,
                                      sdata.column,
                                      rdata.column,
                                      sid.column, 
                                      rid.column,
                                      id.target,   ## for easy function mapping
                                      x.column,
                                      y.column,
                                      radii = 30,
                                      similarity.duration = 168, 
                                      compl.threshold = 90,
                                      direction=c('left','center','right'),
                                      param = 'PM2.5',
                                      new = TRUE) {
  # argument flags
  stopifnot(is.character(date.column), date.column %in% colnames(sdata), date.column %in% colnames(rdata))
  stopifnot(is.character(sdata.column), sdata.column %in% colnames(sdata))
  stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata))
  stopifnot(is.character(sid.column), sid.column %in% colnames(sdata))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  stopifnot(is.character(x.column),x.column %in% colnames(sdata), x.column %in% colnames(rdata))
  stopifnot(is.character(y.column),y.column %in% colnames(sdata), y.column %in% colnames(rdata))
  stopifnot(is.numeric(radii))
  stopifnot(is.numeric(similarity.duration), !is.null(similarity.duration))
  stopifnot(is.numeric(compl.threshold), !is.null(compl.threshold))
  stopifnot(is.numeric(similarity.duration), similarity.duration > 0, !is.null(similarity.duration))
  stopifnot(is.numeric(compl.threshold), compl.threshold > 0, !is.null(compl.threshold))
  stopifnot(direction %in% c('left','center','right'))
  #
  compl.duration <- ifelse(compl.threshold < 1, round(compl.threshold*similarity.duration), compl.threshold)
  radii <- ifelse(is.null(radii), 30, radii)
  radii <- radii*1000
  # target location
  sdata = filter(sdata, !!sym(sid.column) %in% id.target)  ## for easy function mapping
  target.data <- sdata %>% select(!!sym(date.column),!!!syms(sdata.column),
                                  !!sym(x.column),!!sym(y.column))
  
  target.info <- sdata %>% select(!!sym(sid.column),!!sym(x.column),!!sym(y.column)) %>% unique()
  target <- st_as_sf(target.info, crs = 4326L, coords = c(x.column,y.column))
  
  # neighbors reference data
  obs.data <- pull(sdata, sdata.column)
  ind.nona <- which(!is.na(obs.data)) 
  if (length(ind.nona) == 0) {
    output.data <- sdata %>% 
      mutate(flag_similarity = 0)
  } else {
    start.nona <- pull(target.data, date.column)[ind.nona[1]]
    end.nona <- pull(target.data, date.column)[ind.nona[length(ind.nona)]]
    target.data <- target.data %>% 
      filter(!!sym(date.column) >= start.nona & !!sym(date.column) <= end.nona)
    if (sum(!is.na(pull(target.data, sdata.column))) < compl.threshold) {
      output.data <- target.data %>% 
        mutate(flag_similarity = 0)
    } else {
      neighbors.data <- rdata %>% 
        filter(!!sym(date.column) >= start.nona & !!sym(date.column) <= end.nona)
      if (nrow(neighbors.data) == 0) {
        output.data <- sdata %>% 
          mutate(flag_similarity = 1)
      } else {
        neighbors.info <- neighbors.data %>% select(!!sym(rid.column),!!sym(x.column),!!sym(y.column)) %>% unique()
        neighbors <- st_as_sf(neighbors.info, crs = 4326L, coords = c(x.column,y.column))
        distances <- round(st_distance(target,neighbors),1)
        neighbors.info <- neighbors.info %>% 
          mutate(dij = as.vector(distances))
        neighbors.info <- neighbors.info %>% 
          filter(dij >= 0 & dij <= radii)
        #
        if (nrow(neighbors.info) == 0) {
          output.data <- sdata %>% 
            mutate(flag_similarity = 1)
        } else {
          neighbors.data <- neighbors.data %>% 
            filter(!!sym(rid.column) %in% pull(neighbors.info, rid.column))
          neighbors.data.mat <- tidyr::pivot_wider(neighbors.data %>% select(all_of(c(date.column, rid.column, rdata.column))), 
                                                   names_from = !!sym(rid.column), 
                                                   values_from = !!sym(rdata.column)) %>% as_tibble()
          target.neighbors.data.mat <- plyr::join(target.data %>% select(all_of(c(date.column, sdata.column))), 
                                                  neighbors.data.mat, by=date.column) %>% as_tibble()
          neighbors.data.mat <- target.neighbors.data.mat %>% 
            select(-all_of(c(date.column, sdata.column)))
          target.data.mat <- replicate(ncol(neighbors.data.mat), pull(target.neighbors.data.mat, sdata.column))
          #
          # statistical similarity
          similarity.data.mat <- (target.data.mat - neighbors.data.mat)^2
          similarity.stat <- apply(similarity.data.mat, 2, 
                                   function(x) {
                                     similarity <- zoo::rollapply(x, similarity.duration, FUN=function(y) {
                                       sqrt(sum(y, na.rm=TRUE))},by.column= FALSE,fill = NA,align='right')})
          
          neighbors.completeness.stat <- apply(neighbors.data.mat, 2, 
                                               function(x) {
                                                 completeness <- zoo::rollapply(x, similarity.duration, FUN=function(y) {
                                                   sum(!is.na(y))},by.column= FALSE,fill = NA,align='right')})
          
          target.completeness.stat <- zoo::rollapply(pull(target.neighbors.data.mat, sdata.column), similarity.duration, 
                                                     FUN=function(y) {
                                                       sum(!is.na(y))},by.column= FALSE,fill = NA,align='right')
          
          similarity.stat <- ifelse(neighbors.completeness.stat < compl.threshold, NA, similarity.stat)
          neighbors.completeness.stat <- ifelse(neighbors.completeness.stat < compl.threshold, 0, 1)
          #
          target.completeness.flag <- ifelse(target.completeness.stat < compl.threshold, 0, 1)
          neighbors.completeness.flag <- apply(neighbors.completeness.stat, 1, function(x) {
            ifelse(sum(is.na(x)) == length(x), NA, 
                   ifelse(sum(!is.na(x)) >= 1, 1, 0))
          })
          #
          # expected similarity
          expected.similarity.tbl <- target.neighbors.data.mat %>% 
            mutate(Season := quarter(as.POSIXct(!!sym(date.column), "UTC"), fiscal_start = 12)) %>% 
            select(-all_of(c(sdata.column)))
          expected.similarity.tbl <- reshape2::melt(expected.similarity.tbl, id.vars = c(date.column,'Season')) %>% as_tibble() %>% 
            select(-value) %>% 
            rename(!!sym(rid.column) := variable)
          expected.similarity.tbl <- plyr::join(expected.similarity.tbl, neighbors.info %>% select(all_of(c(rid.column, 'dij'))), by=rid.column) %>% as_tibble()
          expected.similarity.tbl <- expected.similarity.tbl %>% 
            mutate(expected_similarity = expected_similarity_fun(season = Season, distance = dij/1000, relax.factor = 0.25))
          expected.similarity.mat <- tidyr::pivot_wider(expected.similarity.tbl %>% select(all_of(c(date.column, rid.column, 'expected_similarity'))), names_from = !!sym(rid.column), values_from = expected_similarity) %>% as_tibble()
          expected.similarity.mat <- expected.similarity.mat %>% 
            select(-!!sym(date.column))
          expected.similarity.stat <- apply(expected.similarity.mat, 2, function(x) {
            similarity <- zoo::rollapply(x, similarity.duration, FUN=function(y) {
              min(DescTools::Mode(y))},by.column= FALSE,fill = NA,align='right')})
          #
          similarity.comparison.mat <- similarity.stat < expected.similarity.stat
          similarity.flag <- apply(similarity.comparison.mat, 1, function(x) {
            flag <- ifelse(sum(is.na(x)) == length(x), NA, sum(x,na.rm=TRUE))
            flag <- ifelse(flag >= 1, 1, 0)
            return(flag)
          })
          #
          # similarity check comparison
          similarity.comparison.mat <- similarity.stat < expected.similarity.stat
          similarity.comparison.mat[is.na(similarity.comparison.mat) | similarity.comparison.mat == FALSE] <- NA
          flag_similarity <- apply(similarity.comparison.mat, 1, function(x) {
            flag <- ifelse(sum(x, na.rm = TRUE) >= 1, 3, 2)
            return(flag)
          })
          flag_similarity <- ifelse(target.completeness.flag == 0, 0, 
                                     ifelse(neighbors.completeness.flag == 0, 1, flag_similarity))
          output.data <- sdata %>% 
            filter(!!sym(date.column) >= start.nona & !!sym(date.column) <= end.nona) %>% 
            mutate(flag_similarity = flag_similarity)
          #
          output.na.startend <- sdata %>% 
            filter(!!sym(date.column) < start.nona |!!sym(date.column) > end.nona) %>% 
            mutate(flag_similarity = NA)
          output.data <- bind_rows(output.na.startend,output.data) %>% 
            as_tibble() 
          output.data <- output.data %>% 
            arrange(!!sym(date.column))
        }
      } 
    }
  }
  output.data <- output.data %>% 
    mutate(flag_similarity = ifelse(is.na(!!sym(sdata.column)),NA,flag_similarity))
  #
  if (new == TRUE) {
    output.data <- output.data %>% 
      mutate(!!sym(paste0(param, '_similarity')) := ifelse(flag_similarity == 0 | flag_similarity == 1 | flag_similarity == 2 | is.na(flag_similarity), NA, obs.data))
  }
  output.data <- output.data %>% 
    as_tibble()
  # 
  return(output.data)
}
#
#
# Correction of sensor data using nearby official measurements
# Development of a calibration function for each day separately using a 30-day rolling window
# Quality flags: 
#   0: No corrected value
#   1: No adequate official sensor -> official measurement site does not exist within  18 km 
#   2: Not adequate data -> data completeness within the temporal window (user defined) is below the completeness threshold (user defined) 
#   3: adequate data -> data completeness within the temporal window (user defined) is above tge completeness threshold (user defined)
calibration_spatial_process <- function(sdata,
                                        rdata,
                                        mdata,
                                        date.column, 
                                        sdata.column,
                                        sdatac.column,
                                        rdata.column,
                                        mdata.column,
                                        sid.column, 
                                        rid.column,
                                        mid.column,
                                        mstation = TRUE,
                                        x.column, 
                                        y.column,
                                        calib.window = 30, 
                                        m.radii = 25,
                                        nmin = 5,
                                        compl.threshold = 100,...) {
  # argument flags
  stopifnot(is.character(date.column), date.column %in% colnames(sdata), date.column %in% colnames(rdata))
  stopifnot(is.character(sdata.column), sdata.column %in% colnames(sdata))
  stopifnot(is.character(sdatac.column), sdatac.column %in% colnames(sdata))
  stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata))
  stopifnot(is.character(mdata.column), mdata.column %in% colnames(mdata))
  stopifnot(is.character(sid.column), sid.column %in% colnames(sdata))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  stopifnot(is.character(mid.column), mid.column %in% colnames(mdata))
  stopifnot(is.logical(mstation))
  stopifnot(is.character(x.column),x.column %in% colnames(sdata), x.column %in% colnames(rdata), x.column %in% colnames(mdata))
  stopifnot(is.character(y.column),y.column %in% colnames(sdata), y.column %in% colnames(rdata), y.column %in% colnames(mdata))
  stopifnot(is.numeric(compl.threshold), !is.null(compl.threshold))
  stopifnot(is.numeric(calib.window),!is.null(calib.window))
  m.radii <- ifelse(is.null(m.radii), 30, m.radii)
  nmin <- ifelse(is.null(nmin), 24, nmin)
  #
  # target location
  target.data <- sdata %>% 
    select(any_of(c(date.column, sdatac.column, sdata.column, mdata.column, sid.column, x.column, y.column)))
  if (mstation) {
    #target.data <- target.data %>% 
    #  mutate(across(any_of(mdata.column), ~NA))
    target.data <- target.data %>% 
      select(-any_of(mdata.column))
  }
  days.intervals <- lubridate::ymd_hms(paste(unique(lubridate::as_date(pull(target.data, !!sym(date.column)))), '23:00:00'))
  target.temporal.windows <- runner::runner(
    x = data.frame(target.data),
    f = function(x) x %>% as_tibble(), 
    k = paste(calib.window,'days'),
    idx = date.column,
    at = days.intervals,
    simplify = FALSE
  )
  days <- lubridate::as_date(as.vector(unlist(lapply(target.temporal.windows, function(x) lubridate::as_date(tail(pull(x,date.column),1))))))
  #
  target.temporal.tbl.list <- vector('list', length(target.temporal.windows))
  nearest.information.tbl.list <- vector('list', length(target.temporal.windows))
  for (i in 1:length(target.temporal.windows)) {
    target.temporal.data <- target.temporal.windows[[i]]
    target.temporal.compl <- target.temporal.data %>% 
      tidyr::drop_na(!!sym(sdata.column)) %>% 
      summarise(n_hours = sum(!is.na(!!sym(sdata.column))))
    r.radii <- lut_distance_fun(season=min(DescTools::Mode(quarter(pull(target.temporal.data, date.column),fiscal_start = 12))))    
    #
    target.temporal.day.tbl <- target.temporal.data %>% 
      mutate(Date = lubridate::as_date(!!sym(date.column))) %>% 
      filter(Date %in% days[i])
    #
    if (pull(target.temporal.compl, 'n_hours') < compl.threshold) {
      target.temporal.day.tbl <- target.temporal.day.tbl %>% 
        mutate(flag_calibration = 0,
                      !!sym(paste0(sdatac.column,'_cor')) := NA,
                      !!sym(paste0(sdatac.column,'_cor_lower')) := NA,
                      !!sym(paste0(sdatac.column,'_cor_upper')) := NA,
                      !!sym(rdata.column) := NA)
    } else {
      # proximal reference station
      target.reference.data <- find_neighbors_fun(sdata=target.temporal.data,
                                                  rdata=rdata,
                                                  date.column = date.column,
                                                  sdata.column = sdata.column,
                                                  rdata.column = rdata.column,
                                                  sid.column = sid.column,
                                                  rid.column = rid.column,
                                                  x.column = x.column, 
                                                  y.column = y.column,
                                                  maxradii = r.radii, 
                                                  maxneigh = NULL,
                                                  compl.threshold = compl.threshold,
                                                  type = 'ref',
                                                  fun='distance')
      nearest.reference <- target.reference.data$neighbors
      target.reference.data <- target.reference.data$data
      #
      # proximal meteorological station
      target.meteo.data <- find_neighbors_fun(sdata=target.temporal.data,
                                              rdata=mdata,
                                              date.column = date.column,
                                              sdata.column = sdata.column,
                                              rdata.column = mdata.column,
                                              sid.column = sid.column,
                                              rid.column = mid.column,
                                              x.column = x.column, 
                                              y.column = y.column,
                                              maxradii = m.radii, 
                                              maxneigh = NULL,
                                              compl.threshold = compl.threshold,
                                              type = 'meteo',
                                              fun = 'distance')
      nearest.meteo <- target.meteo.data$neighbors
      target.meteo.data <- target.meteo.data$data
      #
      if (is.null(nearest.reference)) {
        #
        target.temporal.day.tbl <- target.temporal.day.tbl %>% 
          mutate(flag_calibration = 1,
                        !!sym(paste0(sdatac.column,'_cor')) := NA, 
                        !!sym(paste0(sdatac.column,'_cor_lower')) := NA,
                        !!sym(paste0(sdatac.column,'_cor_upper')) := NA,
                        !!sym(rdata.column) := NA)
      } else {
        # calculation of the minimum reference values
        diurnal.target.reference.data <- target.reference.data %>% 
          mutate(Hour = hour(!!sym(date.column))) %>% 
          group_by(Hour) %>% 
          summarise(variability = var(!!sym(rdata.column),na.rm=TRUE)) %>% 
          arrange(variability)
        min.var <- sort(pull(head(diurnal.target.reference.data, nmin),'Hour'))
        # calculation of the weighted reference values 
        r.distance <- matrix(pull(nearest.reference, 'dij'), nrow = 1)
        r.w <- ((1000*r.radii - r.distance)/(1000*r.radii* r.distance))^2
        for (j in 1:nrow(r.w)) {
          if (sum(is.infinite(r.w[j, ])) > 0) {
            r.w[j, !is.infinite(r.w[j, ])] <- 0
            r.w[j, is.infinite(r.w[j, ])] <- 1
          }
        }
        reference.temporal.mat <- tidyr::pivot_wider(target.reference.data %>% select(all_of(c(date.column, rid.column, rdata.column))), names_from = !!sym(rid.column), values_from = !!sym(rdata.column)) %>% as_tibble()
        reference.temporal.mat <- plyr::join(target.temporal.data %>% select(!!sym(date.column)), reference.temporal.mat, by=date.column) %>% as_tibble()
        reference.temporal.mat <- reference.temporal.mat %>% 
          select(-!!sym(date.column))
        r.w <- matrix(rep(as.vector(r.w),each=nrow(reference.temporal.mat)), nrow=nrow(reference.temporal.mat))
        ind.r.w.mat <- apply(reference.temporal.mat, 2, function(x) ifelse(is.na(x), FALSE, TRUE))
        r.w[!ind.r.w.mat] <- NA
        reference.weighted <-  round(apply(r.w * reference.temporal.mat, 1, sum, na.rm=TRUE)/apply(r.w, 1, sum, na.rm=TRUE),1)
        target.reference.data <- target.temporal.data %>% 
          mutate(Hour = hour(!!sym(date.column))) %>% 
          mutate(!!sym(rdata.column) := reference.weighted)
        #
        if (is.null(nearest.meteo)) {
          target.calibration.data <- target.reference.data 
        } else {
          # calculation of the weighted meteorological values 
          m.distance <- matrix(pull(nearest.meteo, 'dij'), nrow = 1)
          m.w <- ((1000*m.radii - m.distance)/(1000*m.radii* m.distance))^2
          for (j in 1:nrow(m.w)) {
            if (sum(is.infinite(m.w[j, ])) > 0) {
              m.w[j, !is.infinite(m.w[j, ])] <- 0
              m.w[j, is.infinite(m.w[j, ])] <- 1
            }
          }
          target.meteo.tbl.list <- vector('list', length(mdata.column))
          for (j in 1:length(mdata.column)) {
            meteo.temporal.mat <- tidyr::pivot_wider(target.meteo.data %>% select(all_of(c(date.column, mid.column, mdata.column[j]))), names_from = !!sym(mid.column), values_from = !!sym(mdata.column[j])) %>% as_tibble()
            meteo.temporal.mat <- plyr::join(target.temporal.data %>% select(!!sym(date.column)), meteo.temporal.mat, by=date.column) %>% as_tibble()
            target.meteo.tbl <- meteo.temporal.mat %>% 
              select(!!sym(date.column))
            meteo.temporal.mat <- meteo.temporal.mat %>% 
              select(-!!sym(date.column))
            meteo.temporal.mat <- meteo.temporal.mat %>% 
              mutate_all(~imputeTS::na_interpolation(., option = 'linear'))
            mj.w <- matrix(rep(as.vector(m.w),each=nrow(meteo.temporal.mat)), nrow=nrow(meteo.temporal.mat))
            ind.mj.w.mat <- apply(meteo.temporal.mat, 2, function(x) ifelse(is.na(x), FALSE, TRUE))
            mj.w[!ind.mj.w.mat] <- NA
            meteo.weighted <-  round(apply(mj.w * meteo.temporal.mat, 1, sum, na.rm=TRUE)/apply(mj.w, 1, sum, na.rm=TRUE),1)
            if (sum(!is.na(meteo.weighted)) < compl.threshold) {
              meteo.weighted <- NA
            }
            target.meteo.tbl.list[[j]] <- target.meteo.tbl %>%
              mutate(!!sym(mdata.column[j]) := meteo.weighted)
          }
          target.meteo.tbl <- plyr::join_all(target.meteo.tbl.list, by=date.column) %>% as_tibble()
          target.reference.data <- plyr::join(target.reference.data, target.meteo.tbl, by=date.column) %>% as_tibble()
        }
        #
        # check the availability of meteorological data and remove meteorological variables with non-adequate data
        target.reference.compl <- target.reference.data %>% 
          tidyr::drop_na(c(!!sym(sdata.column), !!sym(rdata.column))) %>%
          summarise(across(any_of(c(rdata.column, mdata.column)), ~ sum(!is.na(.))))
        target.meteo.compl <- target.reference.compl %>% 
          select(any_of(mdata.column)) %>%
          mutate(across(any_of(mdata.column), ~if (all(.x > compl.threshold)) .x))
        mdata.common.column <- intersect(names(target.meteo.compl),names(target.reference.data))
        #
        target.calibration.data <- target.reference.data
        target.calibration.data <- target.reference.data %>% 
          mutate(Hour = hour(!!sym(date.column))) %>% 
          filter(Hour %in% min.var) %>% 
          select(-Hour)
        target.calibration.data <- target.calibration.data %>% 
          tidyr::drop_na(!!sym(sdata.column), !!sym(rdata.column))
        #
        if (nrow(target.calibration.data) < compl.threshold) {
          target.temporal.day.tbl <- target.temporal.day.tbl %>% 
            mutate(flag_calibration = 0,
                          !!sym(paste0(sdatac.column,'_cor')) := NA, 
                          !!sym(paste0(sdatac.column,'_cor_lower')) := NA,
                          !!sym(paste0(sdatac.column,'_cor_upper')) := NA,
                          !!sym(rdata.column) := NA)
        } else {
        #
        # correction process
        target.temporal.day.tbl <- plyr::join(target.temporal.day.tbl, target.reference.data %>% select(any_of(c(date.column, rdata.column, mdata.column))), by=date.column) %>% as_tibble()
        target.temporal.day.tbl <- target.temporal.day.tbl %>%
          mutate(across(any_of(mdata.common.column), ~imputeTS::na_interpolation(., option = 'linear')))
        #
        if (length(mdata.common.column) == 0) {
          target.calibration.data <- target.reference.data %>% 
            mutate(Hour = hour(!!sym(date.column))) %>% 
            filter(Hour %in% min.var) %>%
            select(all_of(c(sdata.column, rdata.column)))
          calibration.model <- correction_function(df=target.calibration.data, y= rdata.column, x=sdata.column, params = NULL, conf.level = 0.95, xnew = sdatac.column, newdata =target.temporal.day.tbl)
        } else {
          target.calibration.data <- target.calibration.data %>% 
            mutate(Hour = hour(!!sym(date.column))) %>% 
            filter(Hour %in% min.var) %>%
            select(any_of(c(sdata.column, rdata.column, mdata.common.column)))
          calibration.model <- correction_function(df=target.calibration.data, y= rdata.column, x=sdata.column, params = mdata.common.column, conf.level = 0.95, xnew = sdatac.column, newdata =target.temporal.day.tbl)
        }
        target.calibration.data <- calibration.model$target
        target.temporal.day.tbl <- calibration.model$new
        calibration.performance <- data.frame(RMSDraw = rmsd(pull(target.calibration.data, rdata.column),pull(target.calibration.data, sdata.column)),
                                              RMSDcor = rmsd(pull(target.calibration.data, rdata.column),pull(target.calibration.data, paste0(sdata.column, '_cor')))) %>% as_tibble()
        # calibration  flag
        flag_calibration <- ifelse(sum(is.na(pull(target.temporal.day.tbl, sdatac.column))) == nrow(target.temporal.day.tbl), NA, 
                                   ifelse(is.na(pull(calibration.performance, 'RMSDcor')) & !is.na(pull(calibration.performance, 'RMSDraw')), 2, 
                                          ifelse(pull(calibration.performance, 'RMSDcor') > pull(calibration.performance, 'RMSDraw'), 2, 3)))
        target.temporal.day.tbl <- target.temporal.day.tbl %>% 
          mutate(flag_calibration = flag_calibration)
        }
      }
      #
      target.temporal.day.tbl <- target.temporal.day.tbl %>% 
        select(all_of(c(date.column, sdata.column, sdatac.column, paste0(sdatac.column,'_cor'), paste0(sdatac.column,'_cor_lower'),paste0(sdatac.column,'_cor_upper'), rdata.column, 'flag_calibration')))
    }
    target.temporal.tbl.list[[i]] <- target.temporal.day.tbl
  }
  output.data <- bind_rows(target.temporal.tbl.list) %>% as_tibble()
  output.data <- plyr::join(sdata, output.data %>% select(-any_of(c(sdata.column, sdatac.column, x.column, y.column, sid.column))), by = c(date.column)) %>% as_tibble()
  #
  return(output.data)
}
#
# Validation function
# This function removes the closest reference location
# It corrects the target using the references in the greater neighborhood 
# it validates the corrected target against the excluded reference location
validation_function <- function(sdata,
                                rdata,
                                mdata,
                                date.column, 
                                sdata.column,
                                sdatac.column,
                                rdata.column,
                                mdata.column,
                                sid.column, 
                                rid.column,
                                mid.column,
                                sid.target, 
                                x.column, 
                                y.column,
                                calib.window = 30,
                                calib.radii = 0.5, 
                                m.radii = 30,
                                nmin = 5,
                                compl.threshold = 100,...) {
  #
  # argument flags
  stopifnot(is.character(date.column), date.column %in% colnames(sdata), date.column %in% colnames(rdata))
  stopifnot(is.character(sdata.column), sdata.column %in% colnames(sdata))
  stopifnot(is.character(sdatac.column), sdatac.column %in% colnames(sdata))
  stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata))
  stopifnot(is.character(mdata.column), mdata.column %in% colnames(mdata))
  stopifnot(is.character(sid.column), sid.column %in% colnames(sdata))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  stopifnot(is.character(mid.column), mid.column %in% colnames(mdata))
  stopifnot(is.character(sid.target), !is.null(sid.target))
  stopifnot(is.character(x.column),x.column %in% colnames(sdata), x.column %in% colnames(rdata), x.column %in% colnames(mdata))
  stopifnot(is.character(y.column),y.column %in% colnames(sdata), y.column %in% colnames(rdata), y.column %in% colnames(mdata))
  stopifnot(is.numeric(compl.threshold), !is.null(compl.threshold))
  stopifnot(is.numeric(calib.window),!is.null(calib.window))
  stopifnot(is.numeric(calib.radii), !is.null(calib.radii))
  m.radii <- ifelse(is.null(m.radii), 25, m.radii)
  nmin <- ifelse(is.null(nmin), 24, nmin)
  #
  calib.radii <- calib.radii*1000
  #
  # target location
  target.data <- sdata %>% 
    filter(!!sym(sid.column) %in% sid.target) %>% 
    select(any_of(c(date.column, sdatac.column, sdata.column, mdata.column, sid.column, x.column, y.column)))
  target.info <- sdata %>% filter(!!sym(sid.column) %in% sid.target) %>% select(!!sym(sid.column),!!sym(x.column),!!sym(y.column)) %>% unique()
  target <- st_as_sf(target.info, crs = 4326L, coords = c(x.column,y.column))
  #
  # neighboring locations
  neighbors.data <- rdata %>% 
    select(!!sym(date.column),!!!syms(rdata.column),!!sym(rid.column)) %>% 
    filter(!!sym(date.column) %in% pull(target.data, date.column)) 
  neighbors.info <- rdata %>% 
    filter(!!sym(date.column) %in% pull(target.data, date.column)) %>% 
    select(!!sym(rid.column),!!sym(x.column),!!sym(y.column)) %>% 
    unique()
  if (nrow(neighbors.info) == 0) {
    return(list(data=NULL, stats = NULL))
  } else {
    neighbors <- st_as_sf(neighbors.info, crs = 4326L, coords = c(x.column,y.column))
    # calculation of separation distance and selection of neighbors
    distances <- round(st_distance(target,neighbors),1)
    neighbors.info <- neighbors.info %>% 
      mutate(dij = as.vector(distances))
    neighbors.calib.info <- neighbors.info %>% 
      filter(dij > 0 & dij <= calib.radii)
    if (nrow(neighbors.calib.info) == 0) {
      return(list(data=NULL, stats = NULL))
    } else {
      neighbors.info <- neighbors.info %>% arrange(dij)
      # closest neighbor
      closest.neighbor.info <- neighbors.info %>% 
        head(1)
      closest.neighbor.data <- rdata %>% 
        filter(!!sym(rid.column) %in% pull(closest.neighbor.info, rid.column))
      # remaining neighbors
      other.neighbors.info <- neighbors.info %>% 
        filter(!(!!sym(rid.column) %in% unique(pull(closest.neighbor.info, rid.column))))
      other.neighbors.data <- rdata %>% 
        filter(!!sym(rid.column) %in% pull(other.neighbors.info, rid.column))
      #
      if (nrow(other.neighbors.info) == 0) {
        return(list(data=NULL, stats = NULL))
      } else {
        target.correction <- calibration_spatial_process(sdata = target.data, rdata = other.neighbors.data, mdata = meteo.tbl,
                                                         date.column = date.column, 
                                                         sdata.column = sdata.column, sdatac.column = sdatac.column, rdata.column = rdata.column, mdata.column = mdata.column,
                                                         sid.column = sid.column, rid.column = rid.column, mid.column = mid.column,
                                                         mstation = TRUE, x.column = x.column, y.column = y.column,
                                                         calib.window = 30, m.radii = 30, nmin = 10, compl.threshold = 90)
        #
        target.corrected.data <- target.correction
        calib.compl <- sum(is.na(pull(target.corrected.data, paste0(sdatac.column, '_cor'))))
        if (calib.compl == nrow(target.corrected.data)) {
          return(list(data=NULL, stats = NULL))
        } else {
          target.corrected.data <- target.corrected.data %>% 
            select(all_of(c(date.column, 'flag_calibration', sid.column, sdatac.column, paste0(sdatac.column, '_cor'))))
          validation.data <- plyr::join(target.corrected.data, 
                                        closest.neighbor.data %>% rename(!!sym(paste0(rid.column,'_near')) := !!sym(rid.column)) %>% select(all_of(c(date.column, paste0(rid.column,'_near'), rdata.column))), 
                                        by = date.column) %>% as_tibble()
          valid.compl <-  sum(is.na(pull(validation.data, paste0(sdatac.column, '_cor')) + pull(validation.data, rdata.column)))
          if (valid.compl == nrow(validation.data)) {
            return(list(data=NULL, stats=NULL))
          } else {
            validation.data <- validation.data %>% 
              mutate(!!paste0(sdatac.column,'_cor_flagged') := ifelse(flag_calibration == 0 & !is.na(!!sym(sdatac.column)), !!sym(sdatac.column), !!sym(paste0(sdatac.column,'_cor'))))
            validation.stats <- validation.data %>% 
              group_by(!!sym(sid.column)) %>% 
              tidyr::drop_na(!!sym(rdata.column)) %>%
              summarise(n_raw = sum(!is.na(!!sym(sdatac.column) + !!sym(rdata.column))),
                               mb_raw = md(!!sym(sdatac.column), !!sym(rdata.column), absolute = FALSE), 
                               nmd_raw = norm_md(!!sym(sdatac.column), !!sym(rdata.column), scale=FALSE, absolute = FALSE),
                               mad_raw = md(!!sym(sdatac.column), !!sym(rdata.column), absolute = TRUE), 
                               nmad_raw = norm_md(!!sym(sdatac.column), !!sym(rdata.column), scale=FALSE, absolute = TRUE),
                               rmsd_raw = rmsd(!!sym(sdatac.column), !!sym(rdata.column)), 
                               nrmsd_raw = norm_rmsd(!!sym(sdatac.column), !!sym(rdata.column), scale=FALSE),
                               nsmd_raw = nmsd(!!sym(sdatac.column), !!sym(rdata.column)),
                               mape_raw = mape(!!sym(sdatac.column), !!sym(rdata.column)),
                               r_raw = cor(!!sym(sdatac.column), !!sym(rdata.column), method='pearson', use='pairwise.complete.obs'),
                               #
                               n_cor = sum(!is.na(!!sym(paste0(sdatac.column,'_cor')) + !!sym(rdata.column))),
                               mb_cor = md(!!sym(paste0(sdatac.column,'_cor')), !!sym(rdata.column), absolute = FALSE), 
                               nmd_cor = norm_md(!!sym(paste0(sdatac.column,'_cor')), !!sym(rdata.column), scale=FALSE, absolute = FALSE),
                               mad_cor = md(!!sym(paste0(sdatac.column,'_cor')), !!sym(rdata.column), absolute = TRUE), 
                               nmad_cor = norm_md(!!sym(paste0(sdatac.column,'_cor')), !!sym(rdata.column), scale=FALSE, absolute = TRUE),
                               rmsd_cor = rmsd(!!sym(paste0(sdatac.column,'_cor')), !!sym(rdata.column)), 
                               nrsmd_cor = norm_rmsd(!!sym(paste0(sdatac.column,'_cor')), !!sym(rdata.column), scale=FALSE),
                               nsmd_cor = nmsd(!!sym(paste0(sdatac.column,'_cor')), !!sym(rdata.column)),
                               mape_cor = mape(!!sym(paste0(sdatac.column,'_cor')), !!sym(rdata.column)),
                               r_cor = cor(!!sym(paste0(sdatac.column,'_cor')), !!sym(rdata.column), method='pearson', use='pairwise.complete.obs'),
                               #
                               n_flagged = sum(!is.na(!!sym(paste0(sdatac.column,'_cor_flagged')) + !!sym(rdata.column))),
                               mb_flagged = md(!!sym(paste0(sdatac.column,'_cor_flagged')), !!sym(rdata.column), absolute = FALSE), 
                               nmd_flagged = norm_md(!!sym(paste0(sdatac.column,'_cor_flagged')), !!sym(rdata.column), scale=FALSE, absolute = FALSE),
                               mad_flagged = md(!!sym(paste0(sdatac.column,'_cor_flagged')), !!sym(rdata.column), absolute = TRUE), 
                               nmad_flagged = norm_md(!!sym(paste0(sdatac.column,'_cor_flagged')), !!sym(rdata.column), scale=FALSE, absolute = TRUE),
                               rmsd_flagged = rmsd(!!sym(paste0(sdatac.column,'_cor_flagged')), !!sym(rdata.column)), 
                               nrsmd_flagged = norm_rmsd(!!sym(paste0(sdatac.column,'_cor_flagged')), !!sym(rdata.column), scale=FALSE),
                               nsmd_flagged = nmsd(!!sym(paste0(sdatac.column,'_cor_flagged')), !!sym(rdata.column)),
                               mape_flagged = mape(!!sym(paste0(sdatac.column,'_cor_flagged')), !!sym(rdata.column)),
                               r_flagged = cor(!!sym(paste0(sdatac.column,'_cor_flagged')), !!sym(rdata.column), method='pearson', use='pairwise.complete.obs'))
            validation.stats <- as.data.frame(cbind(Raw = as.vector(unlist(validation.stats %>% select(contains('raw')))), 
                                                    Calibrated = as.vector(unlist(validation.stats %>% select(contains('cor')))),
                                                    Calibrated_flagged = as.vector(unlist(validation.stats %>% select(contains('flagged'))))))
            validation.stats <- validation.stats %>% as_tibble() %>% 
              mutate(Statistic = c('N','MBE', 'nMBE', 'MAE', 'nMAE', 'RMSE', 'nRMSE', 'nMSD', 'MAPE', 'R'))
            validation.stats <- bind_cols(validation.stats, target.info, closest.neighbor.info %>% rename_with(~paste0(.,'_near')))
            return(list(data = validation.data, stats = validation.stats))
          }
        }
      }
    }
  }
}
#
# definition of the overall quality flag
# If the correlation flag is lower than 3 or is missing -> 'other quality'
# If the correlation flag equals to 3 -> 'good quality'
# If the similarity flag is higher than 2 -> 'high quality'
quality_flag_function <- function(data, data.column, 
                                  flag1 = 'flag_correlation', 
                                  flag2 = 'flag_similarity', 
                                  qc.column = NULL) {
  stopifnot(is.character(data.column), data.column %in% colnames(data))
  
  values = pull(data, data.column)
  if (!is.null(qc.column)){
    flag1 = substr(pull(data, qc.column), 4, 4) |> as.numeric()
    flag2 = substr(pull(data, qc.column), 5, 5) |> as.numeric()
  } else {
    stopifnot(is.character(flag1), flag1 %in% colnames(data))
    stopifnot(is.character(flag2) | is.null(flag2) | flag2 != 'flag_similarity', flag2 %in% colnames(data))
    
    flag1 = pull(data, flag1) |> as.numeric()
    flag2 = pull(data, flag2) |> as.numeric()
  }
  
  quality.flag <- rep(NA, nrow(data))
  Good <- which(!is.na(values) & !is.na(flag1) & flag1 == 3)
  Other <- unique(c(which(!is.na(values) & is.na(flag1)),
                    which(!is.na(flag1) & flag1 < 3)))
  N.A <- which(is.na(values))
  quality.flag[Good] <- 2
  quality.flag[Other] <- 3
  quality.flag[N.A] <- NA
  
  if (!is.null(flag2)) {
   High <- which(!is.na(values) & !is.na(flag2) & flag2 == 3)
   quality.flag[High] <- 1
  }
  output.data <- data %>% mutate(quality_flag = quality.flag)
  return(output.data)
} 

qc_flags = function(x){
  # valid range
  range_validity_qctest(data = x, data.column = dc, #new = F,
                        upper.bound = 998, lower.bound = 0.00001) |> 
    
    # constant value
    constant_value_qctest(date.column = dt, data.column = dc, #new = F,
                          persist.duration = 8, compl.duration = 1, min.variation = 0.1, 
                          direction = "center", method = "range") |> 
    # outlier (sensor)
    outlier_qctest(date.column = dt, data.column = dc, #new = F,
                   outlier.duration = 360, compl.duration = 1, outlier.threshold = 10,
                   direction='center')
}


####

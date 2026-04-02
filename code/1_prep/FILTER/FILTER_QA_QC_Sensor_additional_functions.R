#
# expected correlation function
# Pearson's correlation coefficient vs. separation distance for the four calendar seasons
expected_correlation_fun <- function(season, distance, relax.factor = 0.1) {
  Season <- c(1,2,3,4)
  b0 <- c(0.866920787, 0.851493359, 0.764202409, 0.82053746)
  b1 <- c(-0.00265353,-0.00301283,-0.003258162,-0.00285821)
  b2 <- c(5.735157e-06, 6.257868e-06, 7.885344e-06, 6.338316e-06)
  b3 <- c(-5.216934e-09, -5.111353e-09, -7.204441e-09, -5.53291e-09)
  lookup_correlation_distance <- tibble(Season = Season, b0, b1, b2, b3)
  #
  params <- tibble(Season = season, distance = distance)
  params <- plyr::join(params, lookup_correlation_distance, by = 'Season') %>% as_tibble()
  correlation <- pull(params,'b0') + pull(params,'b1')*pull(params,'distance') + pull(params,'b2')*pull(params,'distance')^2 + pull(params,'b3')*pull(params,'distance')^3
  if (!is.null(relax.factor)) {
    correlation <- correlation - relax.factor*correlation
  }  
  return(correlation)
}
#
# expected correlation function
# Pearson's correlation coefficient vs. separation distance for the four calendar seasons
expected_similarity_fun <- function(season, distance, relax.factor = 0.1) {
  Season <- c(1,2,3,4)
  b0 <- c(68.04585477, 58.51618754, 50.71076079, 57.14827933)
  b1 <- c(0.50238293, 0.452549077, 0.286315948, 0.333297642)
  b2 <- c(-0.001295603, -0.001263275, -0.000851074, -0.000914356)
  b3 <- c(1.24202946478091e-06, 1.21303177254065e-06, 8.4149769556485e-07, 8.78846013819774e-07)
  lookup_similarity_distance <- tibble(Season = Season, b0, b1, b2, b3)
  #
  params <- tibble(Season = season, distance = distance)
  params <- plyr::join(params, lookup_similarity_distance, by = 'Season') %>% as_tibble()
  similarity <- pull(params,'b0') + pull(params,'b1')*pull(params,'distance') + pull(params,'b2')*pull(params,'distance')^2 + pull(params,'b3')*pull(params,'distance')^3 
  if (!is.null(relax.factor)) {
    similarity <- similarity - relax.factor*similarity
  }
  return(similarity)
}
#
# look up table for calibration (similarity) distance 
lut_distance_fun <- function(season) {
  Season = c(1,2,3,4)
  Distance = c(11.4847, 12.71141, 19.9706, 17.08543) # in km
  lookup_distance <- tibble(Season = Season, Distance = Distance)
  dist <- pull(lookup_distance %>% filter(Season %in% season), 'Distance')
  return(dist)
}
#
# find proximal sites
# find the n-nearest neighbors based on a statistical metric
# the statistical metrics is the separation distance (Euclidean distance), the correlation coefficient (Pearson's R correlation coefficient) and the distribution similarity (Earth Mover Distance, EMD)
find_neighbors_fun <- function(sdata,
                               rdata,
                               date.column,
                               sdata.column,
                               rdata.column,
                               sid.column,
                               rid.column,
                               x.column, 
                               y.column,
                               maxradii, 
                               maxneigh,
                               compl.threshold,
                               type = c('ref','meteo'),
                               fun = c('distance','correlation','similarity'),...) {
  # argument flags
  stopifnot(is.character(date.column), date.column %in% colnames(sdata), date.column %in% colnames(rdata))
  stopifnot(is.character(sid.column), sid.column %in% colnames(sdata))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  stopifnot(is.character(sdata.column), sdata.column %in% colnames(sdata))
  stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata))
  stopifnot(is.character(x.column),x.column %in% colnames(sdata), x.column %in% colnames(rdata))
  stopifnot(is.character(y.column),y.column %in% colnames(sdata),  y.column %in% colnames(rdata))
  stopifnot(is.numeric(maxradii), !is.null(maxradii))
  stopifnot(is.numeric(compl.threshold), !is.null(compl.threshold))
  stopifnot(type %in% c('ref','meteo'))
  stopifnot(fun %in% c('distance','correlation', 'similarity'))
  #
  # maximum radii in m
  maxradii <- maxradii*1000
  # target location
  target.data <- sdata %>% 
    dplyr::select(!!sym(date.column),!!!syms(sdata.column),!!sym(sid.column),!!sym(x.column),!!sym(y.column))
  target.info <- sdata %>% 
    dplyr::select(!!sym(sid.column),!!sym(x.column),!!sym(y.column)) %>% unique()
  target <- st_as_sf(target.info, crs = 4326L, coords = c(x.column,y.column))
  # neighbors reference data
  neighbors.data <- rdata %>% 
    dplyr::select(!!sym(date.column),!!!syms(rdata.column),!!sym(rid.column)) %>% 
    filter(!!sym(date.column) %in% pull(target.data, date.column)) 
  neighbors.info <- rdata %>% 
    filter(!!sym(date.column) %in% pull(target.data, date.column)) %>% 
    dplyr::select(!!sym(rid.column),!!sym(x.column),!!sym(y.column)) %>% 
    unique()
  if (nrow(neighbors.info) == 0) {
    return(list(data=NULL,neighbors=NULL))
  } else {
    neighbors <- st_as_sf(neighbors.info, crs = 4326L, coords = c(x.column,y.column))
    # calculation of separation distance and selection of neighbors
    distances <- round(st_distance(target,neighbors),1)
    neighbors.info <- neighbors.info %>% 
      mutate(dij = as.vector(distances))
    neighbors.info <- neighbors.info %>% 
      dplyr::arrange(dij) %>% 
      filter(dij > 0 & dij < maxradii)
    if (nrow(neighbors.info) == 0) {
      return(list(data=NULL,neighbors=NULL))
    } else {
      neighbors.info <- neighbors.info %>% dplyr::arrange(dij)
      neighbors.data <- neighbors.data %>% 
        filter(!!sym(rid.column) %in% pull(neighbors.info, rid.column))
      target.neighbors.data <- plyr::join(neighbors.data, target.data %>% dplyr::select(all_of(c(date.column, sdata.column, sid.column))), by=date.column) %>% as_tibble()
      #
      # data completeness and selection of neighboring location 
      target.neighbors.data.compl <- target.neighbors.data %>% 
        dplyr::select(c(!!sym(rid.column), !!!syms(rdata.column),!!sym(sdata.column))) %>% 
        tidyr::drop_na(c(!!!syms(rdata.column),!!sym(sdata.column))) %>%
        group_by(!!sym(rid.column)) %>% 
        dplyr::summarise(across(everything(), ~ sum(!is.na(.x), na.rm = TRUE)))
      target.neighbors.data.compl <- target.neighbors.data.compl %>% 
        filter(!!sym(sdata.column) >= compl.threshold)
      if (nrow(target.neighbors.data.compl) == 0) {
        return(list(data=NULL,neighbors=NULL))
      } else {
        target.neighbors.data.compl <- plyr::join(target.neighbors.data.compl, neighbors.info %>% dplyr::select(!!sym(rid.column),'dij'),by = rid.column) %>% as_tibble() 
        target.neighbors.data.compl <- target.neighbors.data.compl %>% 
          dplyr::arrange(dij)
        # correlation between target and neighboring locations 
        if (type == 'ref') {
          target.neighbors.data.fun <- target.neighbors.data %>% 
            dplyr::select(c(!!sym(rid.column), !!!syms(rdata.column),!!sym(sdata.column))) %>% 
            tidyr::drop_na(c(!!!syms(rdata.column),!!sym(sdata.column))) %>%
            group_by(!!sym(rid.column)) %>% 
            dplyr::summarise(R = cor(!!!syms(rdata.column),!!sym(sdata.column), method = 'pearson', use='pairwise.complete.obs'),
                             W = twosamples::wass_stat(!!!syms(rdata.column),!!sym(sdata.column)))
          #
          target.neighbors.data.compl <- plyr::join(target.neighbors.data.compl, 
                                                    target.neighbors.data.fun, by = rid.column) %>% as_tibble()
          #target.neighbors.data.compl <- target.neighbors.data.compl %>% 
          #  dplyr::arrange(desc(COR))
          #
        }
        if (is.null(maxneigh)) {
          maxneigh <- nrow(target.neighbors.data.compl)
        }
        # selection of the nearest neighboring locations
        if (fun == 'distance') {
          target.neighbors.data <- target.neighbors.data %>% 
            filter(!!sym(rid.column) %in% pull(target.neighbors.data.compl %>% dplyr::arrange(dij) %>% head(maxneigh), rid.column))
        } else if (fun == 'correlation') {
          target.neighbors.data <- target.neighbors.data %>% 
            filter(!!sym(rid.column) %in% pull(target.neighbors.data.compl %>% dplyr::arrange(desc(R)) %>% head(maxneigh), rid.column))
        } else if (fun == 'similarity') {
          target.neighbors.data <- target.neighbors.data %>% 
            filter(!!sym(rid.column) %in% pull(target.neighbors.data.compl %>% dplyr::arrange(W) %>% head(maxneigh), rid.column))
        }
        target.neighbors.info <- neighbors.info %>% 
          filter(!!sym(rid.column) %in% unique(pull(target.neighbors.data,rid.column)))
        #
        return(list(data=target.neighbors.data,neighbors=target.neighbors.info))
      }
    }
  }
}
# 
# correction equation
correction_function <- function(df, y= 'reference', x = 'sensor', params = 'rh', conf.level = 0.9, xnew = NULL, newdata = NULL, robust=TRUE, ...) {
  # argument flags
  stopifnot(is.character(y), y %in% colnames(df))
  stopifnot(is.character(x), y %in% colnames(df))
  stopifnot(is.character(params) | params %in% colnames(df) | is.null(params))
  stopifnot(is.numeric(conf.level), !is.null(conf.level))
  #
  # full formula
  df <- df %>% 
    dplyr::mutate(!!sym(y) := log(!!sym(y) + 0.1),
                  !!sym(x) := log(!!sym(x) + 0.1))
  model.formula <- if (is.null(params)) formula(paste(y, x, sep='~')) else formula(paste(y,paste(x,paste(params, collapse='+'),sep='+'),sep='~'))
  if (robust) {
    initial.model <- MASS::rlm(model.formula, data = data.frame(df), na.action = na.exclude)
  } else {
    initial.model <- lm(model.formula, data = data.frame(df), na.action = na.exclude)
  }
  # remove outliers using MED-MAD function
  df.filtered <- df %>%  
    dplyr::mutate(Bias = residuals(initial.model)) %>%
    filter(abs(Bias - median(Bias, na.rm=TRUE)) < 3*mad(Bias, constant = 1, na.rm=TRUE)) %>% 
    dplyr::select(-Bias)
  #
  # correction model
  if (robust) {
    model <- MASS::rlm(model.formula, data=df.filtered, na.action = na.exclude)
    t_statistic <- coefficients(summary(model))[,3]
    dof <- summary(model)$df[2]
    p.params <- 2*pt(abs(t_statistic), dof, lower.tail = FALSE)
  } else {
    model <- lm(model.formula, data=df.filtered, na.action = na.exclude)
    p.params <- coefficients(summary(model))[,4]
  }
  model.params <- names(p.params)[p.params < (1 - conf.level)]
  model.params <- model.params[!(model.params %in% '(Intercept)')]
  if (length(model.params) == 0) {
    df <- df %>% 
      dplyr::mutate(!!paste0(x,'_cor') := NA,
                    !!paste0(x,'_cor_upper') := NA,
                    !!paste0(x,'_cor_lower') := NA)
    if (!is.null(newdata) | !is.null(xnew)) {
      newdata <- newdata %>% 
        dplyr::mutate(!!paste0(xnew,'_cor') := NA,
                      !!paste0(xnew,'_cor_upper') := NA,
                      !!paste0(xnew,'_cor_lower') := NA) 
    }
  } else {
    model.formula <- formula(paste(y, paste(model.params,collapse='+'), sep='~'))
    if (robust) {
      updated.model <- MASS::rlm(model.formula, data = df.filtered, na.action = na.exclude)
    } else {
      updated.model <- lm(model.formula, data = df.filtered, na.action = na.exclude)
    }
    #
    # predictions
    df.predictions <- predict(updated.model, newdata = df, interval = 'confidence', level = conf.level)
    df.predictions <- df.predictions %>% as_tibble()
    df <- df %>% 
      dplyr::mutate(!!sym(y) := exp(!!sym(y)) - 0.1,
                    !!sym(x) := exp(!!sym(x)) - 0.1) %>%
      dplyr::mutate(!!paste0(x,'_cor') := round(exp(pull(df.predictions, 'fit')) - 0.1, 1),
                    !!paste0(x,'_cor_upper') := round(exp(pull(df.predictions, 'upr')) - 0.1, 1),
                    !!paste0(x,'_cor_lower') := round(exp(pull(df.predictions, 'lwr')) - 0.1, 1))
    if (is.null(newdata) | is.null(xnew)) {
      newdata <- newdata %>% 
        dplyr::mutate(!!paste0(xnew,'_cor') := NA,
                      !!paste0(xnew,'_cor_upper') := NA,
                      !!paste0(xnew,'_cor_lower') := NA) 
    } else {
      newdata <- newdata %>% 
        dplyr::mutate(!!sym(x) := log(!!sym(xnew) + 0.1))
      newdata.predictions <- predict(updated.model, newdata = newdata, interval = 'confidence', level = conf.level)
      newdata.predictions <- newdata.predictions %>% as_tibble() 
      newdata <- newdata %>% 
        dplyr::mutate(!!sym(xnew) := exp(!!sym(x)) - 0.1) %>%
        dplyr::mutate(!!paste0(xnew,'_cor') := round(exp(pull(newdata.predictions, 'fit')) - 0.1, 1),
                      !!paste0(xnew,'_cor_upper') := round(exp(pull(newdata.predictions, 'upr')) - 0.1, 1),
                      !!paste0(xnew,'_cor_lower') := round(exp(pull(newdata.predictions, 'lwr')) - 0.1, 1))
    }
  }
  return(list(target = df, new = newdata))
}
#
# statistical metrics 
# Mean Bias Difference, MBD and Mean Absolute Difference (MAD)
md <- function(x,y, all=FALSE, absolute =FALSE) {
  xy <- data.frame(x,y) 
  if (all == FALSE) {
  xy <- xy %>% drop_na()
  }
  x <- pull(xy, 'x')
  y <- pull(xy, 'y')
  mb <- ifelse(absolute==FALSE, mean(x - y, na.rm=TRUE),mean(abs(x - y), na.rm=TRUE))
  return(mb)
}
# 
# normalized Mean Bias Difference (nMBD)
# the normalization parameters are a) mean and b) standard deviation
norm_md <- function(x, y, scale = FALSE, ...) {
  xy <- data.frame(x,y) 
  xy <- xy %>% drop_na()
  x <- pull(xy, 'x')
  y <- pull(xy, 'y')
  mb <- md(x,y, ...)
  nmb <- ifelse(scale==FALSE, mb/mean(y, na.rm=TRUE), mb/sd(y, na.rm=TRUE))
  return(nmb)
}
# Root Mean Square Deviation (RMSD)
rmsd <- function(x,y, all = FALSE) {
  xy <- data.frame(x,y)
  if (all == FALSE) {
    xy <- xy %>% drop_na()
  }
  x <- pull(xy, 'x')
  y <- pull(xy, 'y')
  rms <- sqrt(mean((x-y)^2,na.rm=TRUE))
  return(rms)
}
# 
# normalized Root Mean Square Difference (nRMSD)
# the normalization parameters are a) mean and b) standard deviation
norm_rmsd <- function(x, y, scale = FALSE, ...) {
  xy <- data.frame(x,y) 
  xy <- xy %>% drop_na()
  x <- pull(xy, 'x')
  y <- pull(xy, 'y')
  rms <- rmsd(x,y, ...)
  nrms <- ifelse(scale==FALSE, rms/mean(y, na.rm=TRUE), rms/sd(y, na.rm=TRUE))
  return(nrms)
}
#
# Mean Absolute Percentafe Error (MAPE)
mape <- function(x, y) {
  xy <- data.frame(x,y) 
  xy <- xy %>% drop_na()
  x <- pull(xy, 'x')
  y <- pull(xy, 'y')
  mpe <- mean(abs(x - y)/abs(y), na.rm=TRUE)
  return(mpe)
}
# 
# normalized standard deviation (nMSD)
nmsd <- function(x, y) {
  xy <- data.frame(x,y) 
  xy <- xy %>% drop_na()
  x <- pull(xy, 'x')
  y <- pull(xy, 'y')
  nms <- (sd(x, na.rm=TRUE) - sd(y, na.rm=TRUE))/sd(y, na.rm=TRUE)
  return(nms)
}

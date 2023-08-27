
#' Generate (#tx) by (#co) distance matrix
#'
#' @param df dataframe
#' @param dist_scaling vector of scaling constants for covariates (NA for exact match)
#' @param metric distance metric to use
#'
#' @return (#tx) by (#co) distance matrix
gen_dm <- function(df,
                   dist_scaling,
                   metric = c("maximum", "euclidean", "manhattan")) {
  metric <- match.arg(metric)

  # NAs indicate exact matches
  #  - hacky workaround: scale exactly-matched covariates by large constant
  #    to stretch any differences to be huge
  # TODO: make this nicer?
  MAX_SCALING <- 1000
  dist_scaling[is.na(dist_scaling)] <- MAX_SCALING

  # pull out df with only covariates
  covs <- df %>% dplyr::select(dplyr::starts_with("X"))
  stopifnot(ncol(covs) == length(dist_scaling))

  # row numbers of each tx/co unit
  tx_obs <- which(as.logical(df$Z))
  co_obs <- which(as.logical(!df$Z))

  # scale covariate df, according to implicit distance metric
  #  - note: coerce non-numeric factor columns to integer values
  #  - note: coerces T/F into 1/2 instead of 1/0, but this is ok
  # covs becomes:  covs * scaling
  #               (nxp)    (pxp)
  #  --> so dist2 is basically using V(x1-x2) rather than (x1-x2),
  #      for V = scaling
  #   - this means that euclidean distance has a (V^T V) scale!
  covs <- covs %>%
    dplyr::mutate(dplyr::across(!dplyr::where(is.numeric),
                                ~as.numeric(as.factor(.)))) %>%
    as.matrix() %*%
    diag(dist_scaling)

  # compute (#tx) x (#co) distance matrix
  # TODO: any way to make this faster?
  dm <- flexclust::dist2(matrix(covs[tx_obs,], ncol = ncol(covs)),   # in case only one tx unit
                         covs[co_obs,],
                         method = metric)
  rownames(dm) <- tx_obs
  colnames(dm) <- co_obs

  return(dm)
}




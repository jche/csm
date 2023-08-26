
#' Generate (#tx) by (#co) distance matrix
#'
#' @param df dataframe
#' @param scaling vector of scaling constants for covariates
#' @param metric distance metric to use
#'
#' @return (#tx) by (#co) distance matrix
gen_dm <- function(df,
                   scaling = 1,
                   metric = c("maximum", "euclidean", "manhattan")) {
  metric <- match.arg(metric)

  # pull out df with only covariates
  covs <- df %>% dplyr::select(dplyr::starts_with("X"))

  # row numbers of each tx/co unit
  tx_obs <- which(as.logical(df$Z))
  co_obs <- which(as.logical(!df$Z))

  # if scaling is scalar c, scale all covariates by c
  if (length(scaling) == 1) {
    scaling <- rep(scaling, ncol(covs))
  }

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
    diag(scaling)

  # compute (#tx) x (#co) distance matrix
  # TODO: any way to make this faster?
  browser()
  # try: dist, rdist, Rfast::Dist
  dm <- flexclust::dist2(matrix(covs[tx_obs,], ncol = ncol(covs)),   # in case only one tx unit
                         covs[co_obs,],
                         method = metric)
  rownames(dm) <- tx_obs
  colnames(dm) <- co_obs

  return(dm)
}

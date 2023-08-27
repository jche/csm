
# functions related to synthetic controls

# main function: generates SC weights
#  - d has tx unit in first row, co in remaining rows
gen_sc_weights <- function(d, dist_scaling,
                           metric = c("maximum", "euclidean", "manhattan")) {
  metric <- match.arg(metric)

  # edge cases
  if (nrow(d) == 0) {
    return(tibble())
  } else if (nrow(d) == 2) {
    return(d %>%
             dplyr::mutate(
               weights = c(1,1)))
  }

  # extreme scales can crash optimizers, so avoid optimizing them
  # TODO: clean this up
  exact_match_cols <- which(is.na(dist_scaling))
  if (length(exact_match_cols) > 0) {
    dist_scaling <- dist_scaling[-exact_match_cols]
  }
  d2 <- d %>%
    dplyr::select(dplyr::starts_with("X")) %>%
    dplyr::select(-dplyr::all_of(exact_match_cols))

  if (metric == "maximum") {
    # run linear program
    sol <- synth_lp(X1 = as.numeric(d2[1,]),
                    X0 = as.matrix(d2[-1,]),
                    V  = diag(dist_scaling))
  } else if (metric == "euclidean") {
    # run osqp
    #  - note: square V, since we use (V^T V) within the euclidean distance!
    sol <- synth_qp(X1 = as.numeric(d2[1,]),
                    X0 = as.matrix(d2[-1,]),
                    V  = diag(dist_scaling^2))
  } else if (metric == "manhattan") {
    stop("Linear program for L1-distance minimization is not currently implemented.")
  }

  d %>%
    mutate(weights = c(1, sol))
}


# Solve the synth QP directly (from ebenmichael/augsynth)
synth_qp <- function(X1, X0, V) {

  Pmat <- X0 %*% V %*% t(X0)
  qvec <- - t(X1) %*% V %*% t(X0)

  n0 <- nrow(X0)
  A <- rbind(rep(1, n0), diag(n0))
  l <- c(1, numeric(n0))
  u <- c(1, rep(1, n0))

  settings = osqp::osqpSettings(verbose = FALSE,
                                eps_rel = 1e-8,
                                eps_abs = 1e-8)
  sol <- osqp::solve_osqp(P = Pmat, q = qvec,
                          A = A, l = l, u = u,
                          pars = settings)

  return(sol$x)
}

# Solve the synth LP
#  - n0 + 1 decision variables: weight for each co unit & one slack variable
synth_lp <- function(X1, X0, V) {

  n0 <- nrow(X0)
  p  <- ncol(X0)

  # build (p x n0) constant matrix
  VX0T <- V %*% t(X0)
  con_mat <- rbind(
    c(0, rep(1, n0)),
    cbind(rep(-1, p), -VX0T),
    cbind(rep(-1, p), VX0T)
  )
  # build (p x 1) constraint vector
  VX1T <- V %*% X1

  obj <- c(1, rep(0,n0))
  dir <- c("=", rep("<=", 2*p))
  rhs <- c(1, -VX1T, VX1T)

  sol <- lpSolve::lp(objective.in = obj,
                     const.mat = con_mat,
                     const.dir = dir,
                     const.rhs = rhs)

  return(sol$solution[-1])
}

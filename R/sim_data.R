
# simple simulated dataset

# generate data mostly in grid (0,1) x (0,1)
gen_df_adv <- function(nc, nt,
                       f0_sd = 0.1,   # homoskedastic noise
                       f0_fun = function(X1, X2) {1},
                       tx_effect_fun = function(X1, X2) {1}
                       # f0_fun = function(X1, X2) { abs(X1-X2) },
                       # tx_effect = function(X1, X2) {(X1-0.5)^2+(X2-0.5)^2},
) {

  SD <- 0.1

  # tx units clustered at (0.25,0.25) and (0.75,0.75)
  dat_txblobs <- dplyr::tibble(
    X1 = c(rnorm(nt/2, mean=0.25, sd=SD),
           rnorm(nt/2, mean=0.75, sd=SD)),
    X2 = c(rnorm(nt/2, mean=0.25, sd=SD),
           rnorm(nt/2, mean=0.75, sd=SD)),
    Z  = T
  )

  # co units clustered at (0.75,0.25) and (0.25,0.75)
  dat_coblobs <- dplyr::tibble(
    X1 = c(rnorm((nc-nt)/2, mean=0.25, sd=SD),
           rnorm((nc-nt)/2, mean=0.75, sd=SD)),
    X2 = c(rnorm((nc-nt)/2, mean=0.75, sd=SD),
           rnorm((nc-nt)/2, mean=0.25, sd=SD)),
    Z  = F
  )

  # some co units uniformly scattered on (0,1) box
  dat_conear <- dplyr::tibble(
    X1 = runif(nt),
    X2 = runif(nt),
    Z  = F
  )

  dat <- dplyr::bind_rows(dat_txblobs, dat_coblobs, dat_conear)

  res <- dat %>%
    dplyr::mutate(
      Y0 = f0_fun(X1,X2) + rnorm(dplyr::n(), mean=0, sd=f0_sd),
      Y1 = f0_fun(X1,X2) + tx_effect_fun(X1,X2) + rnorm(dplyr::n(), mean=0, sd=f0_sd),
      Y  = ifelse(Z, Y1, Y0)) %>%
    dplyr::mutate(id = 1:dplyr::n(), .before=X1)
  # print(paste0("SD of control outcomes: ", round(sd(res$Y0),4)))

  return(res)
}


# generate sample dataset from Hainmueller (2012), exactly
#  - note: generates a big population and samples nt/nco units
# paper settings:
#  - n (total number of units): 300, 600, 1500
#  - r (ratio of co/tx units): 1,2,5
#  - sigma_e: n30 = low overlap, n100 = high overlap, chi5 = weird overlap
#  - outcome: "linear", "nl1", "nl2"
#  - sigma_y: 1
gen_df_hain <- function(nt = 50,
                        nc = 250,
                        sigma_e = c("chi5", "n30", "n100"),
                        outcome = c("linear", "nl1", "nl2"),
                        sigma_y = 1,
                        ATE = 0) {
  sigma_e <- match.arg(sigma_e)
  outcome <- match.arg(outcome)
  NUMSAMP <- max(nt, nc)*10

  if (sigma_e == "n30") {
    eps_e <- rnorm(NUMSAMP, sd=sqrt(30))
  } else if (sigma_e == "n100"){
    eps_e <- rnorm(NUMSAMP, sd=sqrt(100))
  } else if (sigma_e == "chi5") {
    eps_e <- rchisq(NUMSAMP, df=5)
    eps_e <- (eps_e - mean(eps_e)) / sd(eps_e)
    eps_e <- eps_e * sqrt(67.6) + 0.5
  }

  df <- dplyr::as_tibble(mvtnorm::rmvnorm(NUMSAMP,
                          mean = c(0,0,0),
                          sigma = matrix(c(2,1,-1,1,1,-0.5,-1,-0.5,1), ncol=3))) %>%
    dplyr::mutate(V4 = runif(NUMSAMP, min=-3, max=3),
           V5 = rchisq(NUMSAMP, df=1),
           V6 = sample(c(T,F), NUMSAMP, replace=T, prob=c(0.5,0.5)),
           Z  = (V1 + 2*V2 - 2*V3 + V4 - 0.5*V5 + V6 + eps_e) > 0,
           id = 1:dplyr::n()) %>%
    dplyr::relocate(id)

  if (outcome == "linear") {
    df <- df %>%
      dplyr::mutate(Y0 = V1 + V2 + V3 - V4 + V5 + V6 + rnorm(NUMSAMP, sd=sigma_y))
  } else if (outcome == "nl1") {
    df <- df %>%
      dplyr::mutate(Y0 = V1 + V2 + 0.2*V3*V4 - sqrt(V5) + rnorm(NUMSAMP, sd=sigma_y))
  } else if (outcome == "nl2") {
    df <- df %>%
      dplyr::mutate(Y0 = (V1 + V2 + V5)^2 + rnorm(NUMSAMP, sd=sigma_y))
  } else if (outcome == "nl3") {   # not in paper
    df <- df %>%
      dplyr::mutate(Y0 = V1*V2 + V2*V5 + V1*V5 + 3*cos(V3) + 2*sin(V2) +
               V3*V4 + 0.5*V4^2 + rnorm(NUMSAMP, sd=sigma_y))
  }
  df <- df %>%
    dplyr::mutate(Y1 = Y0 + ATE,
           Y  = ifelse(Z, Y1, Y0))

  # get correct number of observations
  if (nt > sum(df$Z) | nc > nrow(df)-sum(df$Z)) {
    warning("Insufficient population pool: increase NUMSAMP")
  }
  df_tx <- df %>%
    dplyr::filter(Z) %>%
    dplyr::slice_sample(n=nt)
  df_co <- df %>%
    dplyr::filter(!Z) %>%
    dplyr::slice_sample(n=nc)

  # bind rows, tx then co
  df_final <- df_tx %>%
    dplyr::bind_rows(df_co)
  names(df_final) <- names(df_final) %>%
    stringr::str_replace("V", "X")

  return(df_final)
}


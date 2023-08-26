
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
  dat_txblobs <- tibble(
    X1 = c(rnorm(nt/2, mean=0.25, sd=SD),
           rnorm(nt/2, mean=0.75, sd=SD)),
    X2 = c(rnorm(nt/2, mean=0.25, sd=SD),
           rnorm(nt/2, mean=0.75, sd=SD)),
    Z  = T
  )

  # co units clustered at (0.75,0.25) and (0.25,0.75)
  dat_coblobs <- tibble(
    X1 = c(rnorm((nc-nt)/2, mean=0.25, sd=SD),
           rnorm((nc-nt)/2, mean=0.75, sd=SD)),
    X2 = c(rnorm((nc-nt)/2, mean=0.75, sd=SD),
           rnorm((nc-nt)/2, mean=0.25, sd=SD)),
    Z  = F
  )

  # some co units uniformly scattered on (0,1) box
  dat_conear <- tibble(
    X1 = runif(nt),
    X2 = runif(nt),
    Z  = F
  )

  dat <- bind_rows(dat_txblobs, dat_coblobs, dat_conear)

  res <- dat %>%
    mutate(Y0 = f0_fun(X1,X2) +
             rnorm(n(), mean=0, sd=f0_sd),
           Y1 = f0_fun(X1,X2) + tx_effect_fun(X1,X2) +
             rnorm(n(), mean=0, sd=f0_sd),
           Y  = ifelse(Z, Y1, Y0)) %>%
    mutate(id = 1:n(), .before=X1)
  # print(paste0("SD of control outcomes: ", round(sd(res$Y0),4)))

  return(res)
}


test_that("matching works", {
  dat <- gen_df_hain(nc = 100, nt = 10)
  res <- csm(dat, dist_scaling = c(1,1,1,1/2,1,NA))
})

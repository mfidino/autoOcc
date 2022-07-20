
all_funcs <- list.files(
  "./R/",
  full.names = TRUE
)

sapply(
  all_funcs,
  source
)

(autologistic_model_fit <- auto_occ(~1 ~1, y=y))
m1 <-  auto_occ(~1 ~1, y=y)
m2 <- auto_occ(~1 ~x, y=y, occ_covs = data.frame(x = rnorm(50)))



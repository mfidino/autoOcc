test_that(
  "vcov.auto_occ_fit",
  {
    data("opossum_det_hist")
    odh <- opossum_det_hist
    opossum_y <- format_y(odh, 1,2,3:6,FALSE)
    opossum_y <- opossum_y[1:30,,]
    opossum_y <- opossum_y[-c(4,12,17),,]
    m1 <- auto_occ(~1~1, y = opossum_y)
    expect_silent(
      vcov(m1)
    )
    expect_true(
      is.matrix(vcov(m1))
    )
    expect_true(
      ncol(vcov(m1)) == nrow(m1@estimates) &
      nrow(vcov(m1) == nrow(m1@estimates))
    )
    expect_silent(
      vcov(m1, type = "psi")
    )
    expect_silent(
      vcov(m1, type = "rho")
    )
    m1@opt$hessian <- NULL
    expect_error(
      vcov(m1, type = "psi")
    )
  }
)

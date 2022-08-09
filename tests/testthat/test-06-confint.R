test_that(
  "confint.auto_occ_fit",
  {
    data("opossum_det_hist")
    odh <- opossum_det_hist
    opossum_y <- format_y(odh, 1,2,3:6,FALSE)
    opossum_y <- opossum_y[1:30,,]
    opossum_y <- opossum_y[-c(4,12,17),,]
    m1 <- auto_occ(~1~1, y = opossum_y)
    expect_silent(
      confint(m1, type = "psi")
    )
    expect_silent(
      confint(m1, type = "rho")
    )
    expect_true(
      is.data.frame(confint(m1, type = "psi"))
    )
    expect_error(
      confint(m1)
    )
    expect_error(
      confint(m1, type = "rho", level = 2)
    )
    expect_error(
      confint(m1, type = "rho", level = "cat")
    )
    expect_error(
      confint(m1, type = "rho", parm = "Impervious")
    )


  }
)

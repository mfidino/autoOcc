test_that(
  "summary",
  {
    data("opossum_det_hist")
    odh <- opossum_det_hist
    opossum_y <- format_y(odh, 1,2,3:6,FALSE)
    opossum_y <- opossum_y[1:30,,]
    opossum_y <- opossum_y[-c(4,12,17),,]
    m1 <- auto_occ(~1~1, y = opossum_y)
    expect_output(
      summary(m1)
    )
    expect_true(
      inherits(
        expect_output(
          summary(m1)
        ),
      "auto_occ_summary"
      )
    )

  }
)

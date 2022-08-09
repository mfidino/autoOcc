test_that(
  "compare_models",
  {
    data("opossum_det_hist")
    data("opossum_covariates")
    odh <- opossum_det_hist
    oc <- opossum_covariates
    # only grabbing the first 30 data.points from each season
    # for this example
    odh <- split(
      odh,
      factor(odh$Season, levels = unique(odh$Season))
    )
    odh <- do.call(
      "rbind.data.frame",
      lapply(
        odh,
        function(x) head(x,30)
      )
    )
    oc <- head(oc,30)

    # function to generate detection history
    opossum_y <- autoOcc::format_y(
      x = odh,
      site_column = "Site",
      time_column = "Season",
      history_columns = "^Week", # start with Week
      report = FALSE
    )


    # scale the covariates (with base R)
    oc_scaled <- as.data.frame(
      lapply(
        oc,
        function(x){
          if(is.numeric(x)){
            scale(x)
          }else{
            x
          }
        }
      )
    )
    # drop site names
    oc_scaled <- oc_scaled[,-1]
    # so many warnings
    m1 <- suppressWarnings(
      auto_occ(
        ~1~1,
        y = opossum_y
      )
    )
    m2 <- suppressWarnings(
      auto_occ(
        ~Impervious~Impervious,
        y = opossum_y,
        det_covs = oc_scaled,
        occ_covs = oc_scaled
      )
    )
    m3 <- suppressWarnings(
      auto_occ(
        ~Impervious + Income~Impervious + Income,
        y = opossum_y,
        det_covs = oc_scaled,
        occ_covs = oc_scaled
      )
    )
    m4 <- suppressWarnings(
      auto_occ(
        ~Income~Income,
        y = opossum_y,
        det_covs = oc_scaled,
        occ_covs = oc_scaled
      )
    )

    expect_silent(
      compare_models(list(m1, m2, m3, m4))
    )
    expect_true(
      all(
        compare_models(list(m1, m2, m3, m4))$model == c("m3", "m4", "m2","m1")
      )
    )
    expect_true(
      all(
        compare_models(
          list(
            intercept = m1,
            impervious = m2,
            impervious_income = m3,
            income = m4
          )
        )$model == c(
          "impervious_income",
          "income",
          "impervious",
          "intercept"
        )
      )
    )
    expect_error(
      compare_models(m1,m2,m3)
    )
    expect_error(
      compare_models(
        list(m1,m2),
        add_formula = "cat"
      )
    )
    expect_error(
      compare_models(
        list(m1,m2),
        digits = "d"
      )
    )
    expect_error(
      compare_models(
        list(m1,m2),
        digits = 0
      )
    )
    expect_error(
      compare_models(
        list(m1,m2),
        digits = -5
      )
    )
    expect_error(
      compare_models(
        list(m1,m1)
      )
    )

    m5 <- suppressWarnings(
      auto_occ(
        ~1~1,
        y = opossum_y[1:20,,]
      )
    )
    expect_error(
      compare_models(
        list(m5, m4, m3)
      )
    )
    tmp_y <- opossum_y
    tmp_y[3,3,3] <- 1
    m6 <- suppressWarnings(
      auto_occ(
        ~1~1,
        y = tmp_y
      )
    )
    expect_error(
      compare_models(
        list(m6, m2,m3)
      )
    )

  }
)

test_that(
  "auto_occ",
  {
    data("opossum_det_hist")
    data("opossum_covariates")
    odh <- opossum_det_hist
    oc <- opossum_covariates

    # function to generate detection history
    opossum_y <- autoOcc::format_y(
      x = odh,
      site_column = "Site",
      time_column = "Season",
      history_columns = "^Week", # start with Week
      report = FALSE
    )
    # removing some data so the model fits faster
    opossum_y <- opossum_y[1:30,,]

    # scale the covariates (with base R)
    oc_scaled <- as.data.frame(
      lapply(
        oc[1:30,],
        function(x){
          if(is.numeric(x)){
            scale(x)
          }else{
            x
          }
        }
      )
    )
    # so many warnings
    expect_warning(
      expect_warning(
        expect_warning(
          m1 <- auto_occ(
            formula = ~Impervious ~ Impervious,
            y = opossum_y,
            det_covs = oc_scaled,
            occ_covs = oc_scaled
          )
       )
     )
    )
    # drop sites with no data, refit model
    new_y <- opossum_y[-c(4,12,17),,]
    new_oc <- oc_scaled[-c(4,12,17),]
    new_oc <- new_oc[,-1]
    expect_silent(
      m1 <- auto_occ(
        formula = ~Impervious ~ Impervious,
        y = new_y,
        det_covs = new_oc,
        occ_covs = new_oc
      )
    )
    # intercept only model
    expect_silent(
      m1 <- auto_occ(
        formula = ~1~1,
        y = new_y,
        det_covs = new_oc,
        occ_covs = new_oc
      )
    )
    # only single formula input
    expect_error(
      m1 <- auto_occ(
        formula = ~1,
        y = new_y,
        det_covs = new_oc,
        occ_covs = new_oc
      )
    )
    # wrong covariate
    expect_error(
      m1 <- auto_occ(
        formula = ~critical~hit,
        y = new_y,
        det_covs = new_oc,
        occ_covs = new_oc
      )
    )
    expect_error(
      m1 <- auto_occ(
        formula = ~1~1,
        y = new_y,
        det_covs = new_oc,
        occ_covs = new_oc,
        level = 1
      )
    )
    expect_error(
      m1 <- auto_occ(
        formula = ~1~1,
        y = new_y,
        det_covs = new_oc,
        occ_covs = new_oc,
        level = "battery"
      )
    )
    expect_silent(
      m1 <- auto_occ(
        formula = ~1~1,
        y = new_y
      )
    )
    expect_silent(
      {
      m1 <- auto_occ(
        formula = ~1~1,
        y = new_y
      )
      tmp_output <- capture.output(
        msum <- summary(m1)
      )
      class(msum) == "summary.auto_occ_fit"
      }
    )
    # just put the whole long detection history in.
    expect_error(
      auto_occ(
        ~1~1,
        y = opossum_det_hist
      )
    )
    expect_error(
      auto_occ(
        ~1~1,
        y = new_y,
        method = "roll initiative"
      )
    )
  }
)

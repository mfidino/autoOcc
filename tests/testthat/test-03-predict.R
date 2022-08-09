test_that(
  "predict",
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
    new_y <- opossum_y[-c(4,12,17),,]
    new_oc <- oc_scaled[-c(4,12,17),]
    new_oc <- new_oc[,-1]
    m1 <- auto_occ(
      ~Impervious ~Impervious + Income,
      y = new_y,
      det_covs = new_oc,
      occ_covs = new_oc
    )
    # site specific predictions rho
    expect_silent(
      predict(m1, type = "rho")
    )
    # same but psi
    expect_silent(
      predict(m1, type = "psi")
    )
    # error because type not included
    expect_error(
      predict(m1)
    )
    ndat <- data.frame(
      Impervious = seq(-3, 3, length.out = 20),
      Income = 0
    )
    # new data predictions psi
    expect_silent(
      predict(m1, type = "psi", newdata = ndat)
    )
    # same but rho
    expect_silent(
      predict(m1, type = "rho", newdata = ndat)
    )
    expect_error(
      predict(m1, type = "psi", newdata = ndat[,1, drop = FALSE])
    )
    expect_error(
      predict(m1, type = "rho", newdata = ndat[,2, drop = FALSE])
    )
    expect_error(
      predict(m1, type = "psi", newdata = ndat, level = "cat")
    )
    expect_error(
      predict(m1, type = "psi", newdata = ndat, level = 1)
    )
    expect_error(
      predict(m1, newdata = ndat)
    )
    expect_error(
      predict(m1, newdata = as.matrix(ndat), type = "psi")
    )
    m2 <- m1
    class(m2) <- "bunny"
    expect_error(
      autoOcc:::predict.auto_occ_fit(m2, newdata = ndat, type = "psi")
    )
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

    season_frame <- list(
      season = matrix(
        factor(odh$Season),
        ncol = dim(opossum_y)[2],
        nrow = dim(opossum_y)[1]
      ),
      Imperv = oc_scaled$Impervious
    )
    m3 <- suppressWarnings(
      auto_occ(
      ~season + Imperv ~season + Imperv,
      y = opossum_y,
      det_covs = season_frame,
      occ_covs = season_frame
    )
    )
    expect_silent(
      predict(m3, type = "psi")
    )
    expect_silent(
      predict(m3, type = "rho")
    )

  }
)

test_that(
  "model_list_estimates",{
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
    mlist <- list(
      suppressWarnings(
        auto_occ(
          ~Impervious ~ Impervious,
          y = opossum_y,
          oc_scaled,
          oc_scaled
        )
      ),
      suppressWarnings(
        auto_occ(
          ~1 ~ 1,
          y = opossum_y,
          oc_scaled,
          oc_scaled
        )
      )
    )

    # should work just fine
    expect_silent(
      model_list_estimates(mlist)
    )
    # retain named list
    names(mlist) <- c("imperv", "null")
    expect_true(
      all(row.names(model_list_estimates(mlist)$ses) == c("imperv","null"))
    )
    # have convergence issues
    bad_y <- opossum_y
    bad_y[,2,] <- 0
    bad_y[,3,] <- 0
    bad_y[,4,] <- 0
    badlist <- mlist
    badlist[[2]]@opt$convergence <- 10
    badlist[[3]] <- mlist[[2]]
    names(badlist)[3] <- "dummy"
    expect_warning(
      model_list_estimates(badlist)
    )

    # vary the data now
    badlist <- mlist
    badlist[[1]]@y[1,1,1] <- 1
    expect_error(
      model_list_estimates(badlist)
    )
    # put in non-list object
    expect_error(
      model_list_estimates(mlist[[1]])
    )
    # put in non auto_occ_Fit
    badlist <- mlist
    badlist[[3]] <- "cranberry"
    expect_error(
      model_list_estimates(badlist)
    )

  }
)

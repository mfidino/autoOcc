test_that(
  "get_dm",
  {
     ###################################
     # Example 1 : Intercept only model
     ###################################

      data("opossum_det_hist")

      # reducing sample size a bit to speed up example run
      opossum_det_hist <- split(
        opossum_det_hist,
        factor(
          opossum_det_hist$Season,
          levels = unique(opossum_det_hist$Season)
         )
       )
      opossum_det_hist <- lapply(
        opossum_det_hist,
        function(x) head(x, 25)
      )
      opossum_det_hist <- do.call(
        "rbind.data.frame",
        opossum_det_hist
      )
      # create y array
      opossum_y <- format_y(
        x = opossum_det_hist,
        site_column = "Site",
        time_column = "Season",
        history_columns = "^Week", # regex for starts with Week
        report = FALSE # defaults to TRUE, turned off for example
      )

      season_frame <- list(
        Season = matrix(
          opossum_det_hist$Season,
          ncol = dim(opossum_y)[2],
          nrow = dim(opossum_y)[1]
        )
      )
      bad_season_frame <- season_frame
      bad_season_frame$Season <- bad_season_frame$Season[,1:3]
      expect_error(
      m1 <- suppressWarnings(
        auto_occ(
         formula =  ~Season  ~Season,
          y = opossum_y,
          det_covs = bad_season_frame,
          occ_covs = bad_season_frame
        )
      )
      )
      expect_error(
        m1 <- suppressWarnings(
          auto_occ(
            formula =  ~Season  ~Season,
            y = opossum_y,
            det_covs = bad_season_frame,
            occ_covs = season_frame
          )
        )
      )
      expect_error(
        m1 <- suppressWarnings(
          auto_occ(
            formula =  ~Season  ~Season + Impervious,
            y = opossum_y,
            det_covs = season_frame,
            occ_covs = season_frame
          )
        )
      )
      expect_error(
        m1 <- suppressWarnings(
          auto_occ(
            formula =  ~Season + Impervious  ~Season,
            y = opossum_y,
            det_covs = season_frame,
            occ_covs = season_frame
          )
        )
      )
       set.seed(3122020)

       # latent occupancy, just doing intercept only.
       nsite <- 25
       nseason <- 6
       nsecondary <- 4

       b0 <- -0.7
       theta <- 1

       occ_prob <- plogis(b0)
       occ_prob_theta <- plogis(b0 + theta)

       z <- matrix(
         NA,
         ncol = nseason,
         nrow = nsite
       )
       z[,1] <- rbinom(
         nsite,
         1,
         occ_prob
       )
       for(t in 2:nseason){
         z[,t] <- rbinom(
           nsite,
           1,
           (1 - z[,t-1]) * occ_prob + z[,t-1] * occ_prob_theta
         )
       }

       # Let's assume we have different observers that go out and
       #  sample the sites, and that precipitation at sites can make
       #  it more difficult to detect the species.



       # Make the precipitation matrix. Each season needs
       # nsecondary columns

       precip <- matrix(
         rnorm(nsite * nseason * nsecondary),
         nrow = nsite,
         ncol = nseason * nsecondary
       )
       # give them more helpful names
       colnames(precip) <- paste0(
         "precip_",
         rep(1:nseason, each = nsecondary),
         "_",
         rep(1:nsecondary, nseason)
       )


       observers <- matrix(
         sample(LETTERS[1:3], nsite * nseason * nsecondary, replace = TRUE),
         nrow = nsite,
         ncol = nseason * nsecondary
       )
       colnames(observers) <- paste0(
         "observer_",
         rep(1:nseason, each = nsecondary),
         "_",
         rep(1:nsecondary, nseason)
       )

       # simulate detection / non-detection data
       y <- array(
         NA,
         dim = c(nsite, nseason, nsecondary)
       )

       # create the regression coefficients
       a0 <- 0.5 #intercept
       a1 <- -0.6 # precip
       a2 <- -1 # observer B
       a3 <- 0.5 # observer C

       alphas <- c(a0,a1,a2,a3)

       for(t in 1:nseason){
         for(j in 1:nsecondary){
           # temporary design matrix
           dmat <- data.frame(
             precip = precip[,grep(paste0("_",t,"_",j), colnames(precip))],
             observer = factor(
               observers[,grep(paste0("_",t,"_",j), colnames(observers))],
               levels = LETTERS[1:3]
             )
           )
           dmat <- model.matrix(~precip + observer, data = dmat)
           logit_prob <- dmat %*% alphas
           my_prob <- plogis(logit_prob)
           y[,t,j] <- rbinom(
             nsite,
             1,
             my_prob * z[,t]
           )
         }
       }
       # make observers a data.frame
       observers <- as.data.frame(observers)

       # make all columns factors with shared levels
       observers <- as.data.frame(
         lapply(
           observers,
           function(x) factor(x, levels = LETTERS[1:3])
         )
       )

       season_frame <- list(
         precip = precip,
         observer = observers
       )

       expect_silent(
       m4 <- auto_occ(
         formula = ~precip + observer ~1,
         y = y,
         det_covs = season_frame
       )
       )


  }
)

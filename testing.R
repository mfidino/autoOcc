# General bookkeeping
nsite <- 25
nyear <- 8
ncovar <- 2
nrepeats <- 5

# for covariates
X <- matrix(
  NA,
  ncol = ncovar,
  nrow = nsite * nyear
)




set.seed(333)
# Create covariates
X <- cbind(1,apply(X, 2, function(x) rnorm(nsite * nyear)))

# Occupancy coefficients, +1 for intercept
psi_betas <- rnorm(ncovar + 1)

# auto-logistic term
theta <- 0.75

# Detection coefficients, decreasing magnitude here
rho_betas <- rnorm(ncovar + 1, 0, 0.5)

# latent state, give same dimensions as X
z <- matrix(NA, ncol = nyear, nrow = nsite)

# Do first year occupancy
psi <- plogis(X %*% psi_betas)
z[,1] <- rbinom(nsite, 1, psi)

# And then the rest, which also uses the theta term
for(year in 2:nyear){
  psi <- plogis(X %*% psi_betas + theta * z[,year-1])
  z[,year] <- rbinom(nsite,1,psi)
}

# Add imperfect detection, make it a matrix with
#  the same dimensions as z. Then multiply by z.
rho <- matrix(
  plogis(X %*% rho_betas),
  ncol = nyear,
  nrow = nsite
) * z

# Create the observed data. Again, same dimensions as z.

y <- array(NA,dim = c(nsite,nyear, nrepeats))
for(i in 1:nrepeats){
  y[,,i] <- rbinom(
    length(rho),
    1,
    rho
  )
}


data_list <- list(
  y=y,
  x = X,
  w = X,
  npar_psi = 4,
  npar_rho = 4,
  nsite = dim(y)[1],
  nseason = dim(y)[2],
  nrep = dim(y)[3]
)


mout <- run.jags("./JAGS/bayes_ao.R",
                 data = data_list,
                 monitor = c("beta","theta","alpha"),
                 n.chains = 4,
                 module = "glm",
                 burnin = 5000,
                 sample = 10000,
                 method = "parallel",
                 inits = list(z = matrix(1, nrow = data_list$nsite, ncol = data_list$nseason)))

oc1 <- data.frame(
  x1 = X[,2],
  x2 = X[,3],
  x3 = X[,4]
)


my_det <- list(
  x1 = oc1$x1,
  x2 = oc1$x2,
  x3 = oc1$x3
)
m2 <- auto_occ(
  formula = ~x1+x2+x3 ~x1+x2+x3 ,
  y=y,
  occ_covs = oc1, det_covs = my_det)

model_summary <- summary(m2)
sl <- round(summary(mout),2)



par(mar = c(3,7,1,1))
plot(1~1, xlim = c(-5,5), ylim = c(0,9), type = "n", bty = "l",
     xlab = "Coefficient estimate",
     ylab = "", yaxt = "n")
axis(2, rev(1:9), model_summary$parameter, las = 2)
for(i in 1:9){
  # add 95% CI and median
  lines(
    x = model_summary[i,c("lower","upper")],
    y = rep(rev(1:9)[i] + 0.15,2),
    col = "blue",
    lwd = 3
  )
  points(
    x = model_summary$Est[i],
    y = rev(1:9)[i] +0.15,
    pch = 21,
    bg = "blue",
    cex = 1.2
  )
  lines(
    x = sl[i,c("Lower95","Upper95")],
    y = rep(rev(1:9)[i] - 0.15,2),
    col = "gray50",
    lwd = 3
  )
  points(
    x = sl[i,"Median"],
    y = rev(1:9)[i] -0.15,
    pch = 21,
    bg = "gray50",
    cex = 1.2
  )
}

 points(y=1:9, x=c(rev(rho_betas), theta, rev(psi_betas)), pch = 21,
        cex = 1.1, bg = "yellow")

legend(
  "topright",
  c("Bayesian model", "Truth", "Frequentist model"),
  pch = 21,
  pt.bg = c("blue","yellow", "gray50"),
  pt.cex = 1.2,
  bty = "n"
)


pdat <- data.frame(x1 = seq(-3,3,0.05), x2=0,x3=0)

model_predict <- predict(m2, newdata = pdat, type = "psi", nsim = 5000)


# do same with bayes model
mcmc <- do.call("rbind", mout$mcmc)[,1:5]



pred1 <- plogis(mcmc %*% t(cbind(1, pdat,0)))
pred2 <- plogis(mcmc %*% t(cbind(1, pdat,1)))

my_est <- pred1 / (pred1 + (1 - pred2))
mcmc_est <- t(apply(my_est, 2, quantile, probs= c(0.025,0.5,0.975)))[,c(2,1,3)]

{
  bbplot::blank(ylim = c(0,1), xlim = c(-3,3), bty = "l", xaxs = "i", yaxs = "i")
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 0.5)
  bbplot::axis_text(side = 2, line = 0.5, las = 1)
  bbplot::axis_text("Covariate", side = 1, line = 2.25, cex = 1.2)
  bbplot::axis_text("Occupancy", side = 2, line = 2.25, cex = 1.2)
  bbplot::ribbon(
    x = pdat$x1,
    y = model_predict[,c("lower","upper")],
    col = "purple",
    alpha = 0.4
  )
  bbplot::ribbon(
    x = pdat$x1,
    y = mcmc_est[,2:3],
    col = "green",
    alpha = 0.4
  )
  lines(x = pdat$x1, y = model_predict$estimate, col = "purple", lwd = 2)

  lines(x = pdat$x1, y = mcmc_est[,1], col = "green", lwd =2, lty = 2)
  legend(
    "topright",
    c("Frequentist", "Bayesian"),
    lty = c(1,2),
    lwd = 4,
    col = c("purple", "green"),
    bty = "n",
    title = "Model type fitted to data")
}


mcmc <- do.call("rbind", mout$mcmc)[,1:5]

pmcmc <- hm[,1:5]

pred1 <- plogis(pmcmc %*% t(cbind(1, pdat,0)))
pred2 <- plogis(pmcmc %*% t(cbind(1, pdat,1)))

my_est <- pred1 / (pred1 + (1 - pred2))
bootsie <- t(apply(my_est, 2, quantile, probs= c(0.025,0.5,0.975)))

range1 <- hm$upper - hm$lower
range2 <- my_est[,3] - my_est[,1]
range3 <- bootsie[,3] - bootsie[,1]

plot(range2 ~ pdat$x1, type = 'l')
lines(range1 ~ pdat$x1, lty = 2)
longshot <- t(apply(pred2, 2, quantile, probs= c(0.025,0.5,0.975)))

plot(mcmc_est[,2] ~ pdat$x1, type = "l", ylim = c(0,1))

lines(bootsie[,2] ~ pdat$x1, col = "red", lwd = 2)

plot(bootsie[,2] ~ pdat$x1, type = "l", ylim = c(0,1))
lines(bootsie[,1] ~ pdat$x1, lty = 2, lwd = 2)
lines(bootsie[,3] ~ pdat$x1, lty = 2, lwd = 2)

plot(mcmc_est[,2] ~ pdat$x1, type = "l", ylim = c(0,1))
lines(mcmc_est[,1] ~ pdat$x1, lty = 2, lwd = 2)
lines(mcmc_est[,3] ~ pdat$x1, lty = 2, lwd = 2)

plot(model_predict$estimate ~ pdat$x1, type = "l", ylim = c(0,1))
lines(model_predict$lower ~ pdat$x1, lty = 2, lwd = 2)
lines(model_predict$upper ~ pdat$x1, lty = 2, lwd = 2)


# use the data examples
data("opossum_det_hist")
data("opossum_covariates")
odh <- opossum_det_hist
oc <- opossum_covariates

# function to generate detection history
opossum_y <- autoOcc::format_y(
  x = odh,
  site_column = "Site",
  time_column = "Season",
  history_columns = "^Week" # start with Week
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
# we are using the same covariates for detection and occupancy
#  so we can use the same data.frame.
m1 <- autoOcc::auto_occ(
  ~Impervious + Income ~ Impervious + Income,
  y = opossum_y,
  det_covs = oc_scaled,
  occ_covs = oc_scaled
)

m2 <- autoOcc::auto_occ(
  ~Impervious  ~ Impervious,
  y = opossum_y,
  det_covs = oc_scaled,
  occ_covs = oc_scaled
)

m3 <- autoOcc::auto_occ(
  ~Income  ~ Income,
  y = opossum_y,
  det_covs = oc_scaled,
  occ_covs = oc_scaled
)

# intercept only model
m4 <- autoOcc::auto_occ(
  ~1~1,
  y = opossum_y
)

my_results <- autoOcc::compare_models(
  list(m1,m2,m3,m4),
  TRUE,
  digits = 2
)

# global model is the best, look at change in occupancy.
# first make the prediction data.frame with a realistic
# range based on the actual data and not the scaled data.
imperv_real <- data.frame(
  Impervious = seq(20,80,0.5),
  Income = 0
)

# then create the scaled data.frame and use that for predictions
imperv <- imperv_real
imperv$Impervious <- (
  imperv$Impervious - mean(oc$Impervious)
  ) / sd(oc$Impervious)

# the model prediction
opo_imperv <- autoOcc::predict(m1, "psi", imperv)

windows(9,4)
par(mfrow = c(1,2))
plot(
  opo_imperv$estimate ~ imperv_real$Impervious,
  bty = "l",
  type = "l",
  las = 1,
  ylab = "Occupancy",
  xlab= "Impervious Cover",
  ylim = c(0,1),
  lwd = 3
)
lines(opo_imperv$lower ~ imperv_real$Impervious, lwd = 2, lty = 2)
lines(opo_imperv$upper ~ imperv_real$Impervious, lwd = 2, lty = 2)


# do the same thing with income
income_real <- data.frame(
  Impervious = 0,
  Income = seq(40000, 160000, by = 500)
)

income <- income_real
income$Income <- (income$Income - mean(oc$Income)) / sd(oc$Income)

opo_income <- autoOcc::predict(m1, "psi", income)

plot(
  opo_income$estimate ~ income_real$Income,
  bty = "l",
  type = "l",
  las = 1,
  ylab = "Occupancy",
  xlab= "Per Capita Income (US Dollar)",
  ylim = c(0,1),
  lwd = 3
)
lines(opo_income$lower ~ income_real$Income, lwd = 2, lty = 2)
lines(opo_income$upper ~ income_real$Income, lwd = 2, lty = 2)



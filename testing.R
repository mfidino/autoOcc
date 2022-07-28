
all_funcs <- list.files(
  "./R/",
  full.names = TRUE
)

sapply(
  all_funcs,
  source
)


{

nsite <- 500
nseason <- 6
nrep <- 4

b0 <- 0.2
b1 <- -0
theta <- 0.7

a0 <- 0.75
a1 <- 0

x <- rnorm(nsite)
z <- matrix(NA, nrow = nsite, ncol = nseason)
y <- array(NA, dim = c(nsite,nseason,nrep))

tmp_psi <- cbind(1,x) %*% c(b0,b1)
tmp_psi2 <- tmp_psi + theta
psi_init <- plogis(tmp_psi)
#psi_init <- plogis(tmp_psi)/(plogis(tmp_psi) + (1 - plogis(tmp_psi2)))

z[,1] <- rbinom(
  nsite,
  1,
  psi_init
)

for(i in 2:nseason){
  psi_year <- plogis(
    cbind(1,x,z[,i-1]) %*% c(b0,b1,theta)
  )
  z[,i] <- rbinom(
    nsite,
    1,
    psi_year
  )
}

for(i in 1:nseason){
  rho_year <- plogis(
    cbind(1,x) %*% c(a0,a1)
  )
  rho_year <- rho_year * z[,i]

  for(j in 1:nsite){
  y[j,i,] <-  rbinom(nrep,1,rho_year[j])
  }
}


my_occ <- data.frame(
  x = x
)

my_det <- list(
  x = matrix(x, nrow = nsite, ncol = nseason)
)
}
m1 <- auto_occ(~1 ~1, y=y)
b0
theta
m1


m3 <- auto_occ(~1 ~x , y=y,
               occ_covs = my_occ, det_covs = my_det)
m4 <- auto_occ(~x ~1, occ_covs = my_occ, det_covs = my_det, y=y)


my_fit <- compare_models(
  list(m1, m2, m3, m4),
  #list(m1 = m1, m2 = m2, m3= m3, m4=m4),
  add_formula = TRUE, digits = 2)

my_fit

nd <- data.frame(x = seq(-3,3,length.out = 200))

my_pred <- predict(m2, newdata = nd, type = "psi")

plot(my_pred$estimate ~ nd$x, type = "l", ylab = "Pr(occupancy)",
     xlab = "Covariate", bty = "l", las = 1, ylim = c(0,1))
lines(my_pred$lower ~ nd$x, lty = 2)
lines(my_pred$upper ~ nd$x, lty = 2)


library(runjags)


round(summary(mout),2)



# General bookkeeping
nsite <- 25
nyear <- 8
ncovar <- 3
nrepeats <- 5

# for covariates
X <- matrix(
  NA,
  ncol = ncovar,
  nrow = nsite
)

set.seed(333)
# Create covariates
X <- cbind(1,apply(X, 2, function(x) rnorm(nsite)))

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
y <- matrix(
  rbinom(
    length(rho),
    nrepeats,
    rho
  ),
  ncol = nyear,
  nrow = nsite
)


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


mout <- run.jags("./R/bayes_ao.R",
                 data = data_list,
                 monitor = c("beta","theta","alpha"),
                 n.chains = 2,
                 module = "glm",
                 method = "parallel",
                 inits = list(z = matrix(1, nrow = data_list$nsite, ncol = data_list$nseason)))

oc1 <- data.frame(
  x1 = X[,2],
  x2 = X[,3],
  x3 = X[,4]
)

my_det <- list(
  x1 = matrix(oc1$x1, nrow = nsite, ncol = nseason),
  x2 = matrix(oc1$x2, nrow = nsite, ncol = nseason),
  x3 = matrix(oc1$x3, nrow = nsite, ncol = nseason)
)
m2 <- auto_occ(
  ~x1+x2+x3 ~x1+x2+x3 ,
  y=y,
  occ_covs = oc1, det_covs = my_det)
jj <- summary(m2)
sl <- round(summary(mout),2)



par(mar = c(3,7,1,1))
plot(1~1, xlim = c(-5,5), ylim = c(0,9), type = "n", bty = "l",
     xlab = "Coefficient estimate",
     ylab = "", yaxt = "n")
axis(2, rev(1:9), jj$parameter, las = 2)
for(i in 1:9){
  # add 95% CI and median
  lines(
    x = jj[i,c("lower","upper")],
    y = rep(rev(1:9)[i] + 0.15,2),
    col = "blue",
    lwd = 3
  )
  points(
    x = jj$Est[i],
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

hm <- predict(m2, newdata = pdat, type = "psi")

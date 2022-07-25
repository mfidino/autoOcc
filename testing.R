
all_funcs <- list.files(
  "./R/",
  full.names = TRUE
)

sapply(
  all_funcs,
  source
)


{set.seed(11911)
nsite <- 50
nseason <- 4
nrep <- 4

b0 <- 0.5
b1 <- -1
theta <- 0.7

a0 <- 0.75
a1 <- -0.3

x <- rnorm(nsite)
z <- matrix(NA, nrow = nsite, ncol = nseason)
y <- array(NA, dim = c(nsite,nseason,nrep))

tmp_psi <- cbind(1,x) %*% c(b0,b1)
tmp_psi2 <- tmp_psi + theta

psi_init <- plogis(tmp_psi)/(plogis(tmp_psi) + (1 - plogis(tmp_psi2)))

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
  my_vals <-  rbinom(nsite*nrep,1,rep(rho_year, each = nrep))
  tmp <- t(y[,i,])
  tmp[] <- my_vals
  tmp <- t(tmp)
  y[,i,] <- tmp
}


my_occ <- data.frame(
  x = x
)

my_det <- list(
  x = matrix(x, nrow = nsite, ncol = nseason)
)
}
m1 <-  auto_occ(~1 ~1, y=y)
m2 <- auto_occ(~x ~x , y=y,
               occ_covs = my_occ, det_covs = my_det)
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

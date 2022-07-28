model{
  for(i in 1:nsite){
    logit(psi[i,1]) <- inprod(
      x[i,],
      beta
    )
    z[i,1] ~ dbern(psi[i,1])
    for(t in 2:nseason){
      logit(psi[i,t]) <- inprod(
        x[i,],
        beta
      ) + theta * z[i,t-1]
      z[i,t] ~ dbern(psi[i,t])
    }
  }
  for(i in 1:nsite){
    for(tt in 1:nseason){
      for(j in 1:nrep){
        logit(rho[i,tt,j]) <- inprod(
          w[i,],
          alpha
        )
        y[i,tt,j] ~ dbern(rho[i,tt,j] * z[i,tt])
      }
    }
  }
    for(ipsi in 1:npar_psi){
      beta[ipsi] ~ dlogis(0,1)
    }
    theta ~ dlogis(0,1)
    for(irho in 1:npar_rho){
      alpha[irho] ~ dlogis(0,1)
    }

}

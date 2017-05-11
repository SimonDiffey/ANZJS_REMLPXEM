#------------------
# Old specification
#------------------
remlpxem_old <- function(y=y, X=X, Z=Z, dd.start=1, sigma2_e.start=1, maxiter=1500, pxem=TRUE)
{
  dd <- lambda <- sigma2_a <- sigma2_e <- NULL
  # Starting values
  dd[1] <- dd.start
  sigma2_e[1] <- sigma2_e.start
  lambda[1] <- 1
  sigma2_a[1] <- lambda[1]^2*dd[1]
  Kappa <- list()
  Kappa.change <- LL <- NULL
  Kappa[[1]] <- c(sigma2_a[1], sigma2_e[1])
  # Form constants
  nn <- dim(X)[1]
  tt <- dim(X)[2]
  bb <- dim(Z)[2]
  K <- diag(nn) - X%*%solve(t(X)%*%X)%*%t(X)
  W <- cbind(X,Z)
  for(w in 1:maxiter)
  { 
    G <- sigma2_a[w]*diag(bb)
    Ginv <- (1/sigma2_a[w])*diag(bb)
    H <- Z%*%G%*%t(Z) + sigma2_e[w]*diag(nn)
    Hinv <- solve(H)
    beta <- solve(t(X)%*%Hinv%*%X)%*%t(X)%*%Hinv%*%y
    P <- Hinv - Hinv%*%X%*%solve(t(X)%*%Hinv%*%X)%*%t(X)%*%Hinv
    S <- (1/sigma2_e[w])*K
    Rinv <- (1/sigma2_e[w])*diag(nn)
    C.ZZ <- solve(t(Z)%*%S%*%Z + Ginv)
    C.XZ <- -solve(t(X)%*%Hinv%*%X)%*%t(X)%*%Hinv%*%Z%*%G
    C_ZZ <- t(Z)%*%Rinv%*%Z + Ginv
    C_XX <- t(X)%*%Rinv%*%X
    C_XZ <- t(X)%*%Rinv%*%Z
    C <- rbind(cbind(C_XX,C_XZ), cbind(t(C_XZ),C_ZZ))
    Cinv <- solve(C)
    u.tilde <- G%*%t(Z)%*%P%*%y
    e.tilde <- sigma2_e[w]*P%*%y
    LL[w] <- -0.5 * ( determinant(t(X)%*%Hinv%*%X, logarithm=TRUE)$modulus + determinant(H, logarithm=TRUE)$modulus + t(y)%*%P%*%y )
    # Updates
    sigma2_e[w+1] <- 1/nn * ( t(e.tilde)%*%e.tilde + sum(diag( W%*%Cinv%*%t(W))) )
    dd[w+1] <- (1/bb)*(t(u.tilde)%*%u.tilde + sum(diag(C.ZZ))) 
    if(pxem==TRUE)
    {
      lambda[w+1] <- (t(u.tilde)%*%t(Z)%*%(y-X%*%beta)-sum(diag(t(Z)%*%X%*%C.XZ)) ) / ( t(u.tilde)%*%t(Z)%*%Z%*%u.tilde + sum(diag(t(Z)%*%Z%*%C.ZZ)) )
    }
    if(pxem==FALSE)
    {
      lambda[w+1] <- 1
    }
    sigma2_a[w+1] <- lambda[w+1]^2*dd[w+1]
    Kappa[[w+1]] <- c(sigma2_a[w+1], sigma2_e[w+1])
    Kappa.change[w] <- sqrt( (t(Kappa[[w+1]] - Kappa[[w]])%*%(Kappa[[w+1]] - Kappa[[w]])) / (t(Kappa[[w]])%*%Kappa[[w]]) )
    if(Kappa.change[w] < 10^-8)
    {
      cat("Converged in",w,"iterations")
      break
    }
  }
  list(Kappa, LL)  
}

#------------------
# New specification
#------------------
remlpxem_new <- function(y=y, X=X, Z=Z, dd.start=1, sigma2_e.start=1, maxiter=1500, pxem=TRUE)
{
  dd <- lambda <- sigma2_a <- sigma2_e <- NULL
  # Starting values
  dd[1] <- dd.start
  sigma2_e[1] <- sigma2_e.start
  lambda[1] <- 1
  sigma2_a[1] <- lambda[1]^2*dd[1]
  Kappa <- list()
  Kappa.change <- LL <- NULL
  Kappa[[1]] <- c(sigma2_a[1], sigma2_e[1])
  # Form constants
  nn <- dim(X)[1]
  tt <- dim(X)[2]
  bb <- dim(Z)[2]
  K <- diag(nn) - X%*%solve(t(X)%*%X)%*%t(X)
  W <- cbind(X,Z)
  for(w in 1:maxiter)
  { 
    G <- sigma2_a[w]*diag(bb)
    Ginv <- (1/sigma2_a[w])*diag(bb)
    H <- Z%*%G%*%t(Z) + sigma2_e[w]*diag(nn)
    Hinv <- solve(H)
    beta <- solve(t(X)%*%Hinv%*%X)%*%t(X)%*%Hinv%*%y
    P <- Hinv - Hinv%*%X%*%solve(t(X)%*%Hinv%*%X)%*%t(X)%*%Hinv
    S <- (1/sigma2_e[w])*K
    Rinv <- (1/sigma2_e[w])*diag(nn)
    C.ZZ <- solve(t(Z)%*%S%*%Z + Ginv)
    C.XZ <- -solve(t(X)%*%Hinv%*%X)%*%t(X)%*%Hinv%*%Z%*%G
    C_ZZ <- t(Z)%*%Rinv%*%Z + Ginv
    C_XX <- t(X)%*%Rinv%*%X
    C_XZ <- t(X)%*%Rinv%*%Z
    C <- rbind(cbind(C_XX,C_XZ), cbind(t(C_XZ),C_ZZ))
    Cinv <- solve(C)
    u.tilde <- G%*%t(Z)%*%P%*%y
    e.tilde <- sigma2_e[w]*P%*%y
    LL[w] <- -0.5 * ( determinant(t(X)%*%Hinv%*%X, logarithm=TRUE)$modulus + determinant(H, logarithm=TRUE)$modulus + t(y)%*%P%*%y )
    # Updates
    sigma2_e[w+1] <- 1/(nn-tt) * ( t(y-Z%*%u.tilde)%*%K%*%(y-Z%*%u.tilde) + sum(diag(t(Z)%*%K%*%Z%*%C.ZZ)) )
    dd[w+1] <- (1/bb)*(t(u.tilde)%*%u.tilde + sum(diag(C.ZZ))) 
    if(pxem==TRUE)
    {
      lambda[w+1] <- t(y)%*%K%*%Z%*%u.tilde / ( t(u.tilde)%*%t(Z)%*%K%*%Z%*%u.tilde + sum(diag(t(Z)%*%K%*%Z%*%C.ZZ)) )
      # # Below is reformulation based on equation 16 in ANZJS paper
      # lambdaP <- (t(u.tilde)%*%t(Z)%*%(y-X%*%beta)-sum(diag(t(Z)%*%X%*%C.XZ)) ) / ( t(u.tilde)%*%t(Z)%*%Z%*%u.tilde + sum(diag(t(Z)%*%Z%*%C.ZZ)) )
      # bbb <- t(u.tilde)%*%t(Z)%*%K%*%Z%*%u.tilde + sum(diag(t(Z)%*%K%*%Z%*%C.ZZ))
      # aaa <- t(u.tilde)%*%t(Z)%*%Z%*%u.tilde + sum(diag(t(Z)%*%Z%*%C.ZZ))
      # lambda[w+1] <- (aaa/bbb)*(lambdaP-1) + 1
    }
    if(pxem==FALSE)
    {
      lambda[w+1] <- 1
    }
    sigma2_a[w+1] <- lambda[w+1]^2*dd[w+1]
    Kappa[[w+1]] <- c(sigma2_a[w+1], sigma2_e[w+1])
    Kappa.change[w] <- sqrt( (t(Kappa[[w+1]] - Kappa[[w]])%*%(Kappa[[w+1]] - Kappa[[w]])) / (t(Kappa[[w]])%*%Kappa[[w]]) )
    if(Kappa.change[w] < 10^-8)
    {
      cat("Converged in",w,"iterations")
      break
    }
  }
  list(Kappa, LL)  
}

#-----------------------------------
# Read data and form design matrices
#-----------------------------------
library(MASS)
library(agridat)

# Lamb weight data
lamb.df <- read.csv("LambWeightData.csv")
lamb.df$Sire <- factor(lamb.df$Sire)
lamb.df$DAMage <- factor(lamb.df$DAMage)
lamb.df$Line <- factor(lamb.df$Line)
str(lamb.df)

yl <- matrix(data=lamb.df$weight, nrow=62)  
Xl <- model.matrix(~1+DAMage+Line, data=lamb.df)
Zl <- model.matrix(~-1+Sire, data=lamb.df)

# Soybean variety trial yield data
str(weiss.incblock)

yw <- matrix(data=weiss.incblock$yield, nrow=186)  
Xw <- model.matrix(~1+gen, data=weiss.incblock)
Zw <- model.matrix(~-1+block, data=weiss.incblock)

#--------------------------------------------------
# Iterations to convergence: Table 2 of ANZJS paper
#--------------------------------------------------

# SV1
remlem.old.lamb.sv1 <- remlpxem_old(y=yl, X=Xl, Z=Zl, dd.start=2, sigma2_e.start=2, maxiter=500, pxem=FALSE)
remlem.new.lamb.sv1 <- remlpxem_new(y=yl, X=Xl, Z=Zl, dd.start=2, sigma2_e.start=2, maxiter=500, pxem=FALSE)
remlpxem.old.lamb.sv1 <- remlpxem_old(y=yl, X=Xl, Z=Zl, dd.start=2, sigma2_e.start=2, maxiter=500, pxem=TRUE)
remlpxem.new.lamb.sv1 <- remlpxem_new(y=yl, X=Xl, Z=Zl, dd.start=2, sigma2_e.start=2, maxiter=500, pxem=TRUE)
remlem.old.soy.sv1 <- remlpxem_old(y=yw, X=Xw, Z=Zw, dd.start=1, sigma2_e.start=1, maxiter=50, pxem=FALSE)
remlem.new.soy.sv1 <- remlpxem_new(y=yw, X=Xw, Z=Zw, dd.start=1, sigma2_e.start=1, maxiter=50, pxem=FALSE)
remlpxem.old.soy.sv1 <- remlpxem_old(y=yw, X=Xw, Z=Zw, dd.start=1, sigma2_e.start=1, maxiter=50, pxem=TRUE)
remlpxem.new.soy.sv1 <- remlpxem_new(y=yw, X=Xw, Z=Zw, dd.start=1, sigma2_e.start=1, maxiter=50, pxem=TRUE)

# SV2
remlem.old.lamb.sv2 <- remlpxem_old(y=yl, X=Xl, Z=Zl, dd.start=3, sigma2_e.start=2, maxiter=500, pxem=FALSE)
remlem.new.lamb.sv2 <- remlpxem_new(y=yl, X=Xl, Z=Zl, dd.start=3, sigma2_e.start=2, maxiter=500, pxem=FALSE)
remlpxem.old.lamb.sv2 <- remlpxem_old(y=yl, X=Xl, Z=Zl, dd.start=3, sigma2_e.start=2, maxiter=500, pxem=TRUE)
remlpxem.new.lamb.sv2 <- remlpxem_new(y=yl, X=Xl, Z=Zl, dd.start=3, sigma2_e.start=2, maxiter=500, pxem=TRUE)
remlem.old.soy.sv2 <- remlpxem_old(y=yw, X=Xw, Z=Zw, dd.start=4, sigma2_e.start=8, maxiter=50, pxem=FALSE)
remlem.new.soy.sv2 <- remlpxem_new(y=yw, X=Xw, Z=Zw, dd.start=4, sigma2_e.start=8, maxiter=50, pxem=FALSE)
remlpxem.old.soy.sv2 <- remlpxem_old(y=yw, X=Xw, Z=Zw, dd.start=4, sigma2_e.start=8, maxiter=50, pxem=TRUE)
remlpxem.new.soy.sv2 <- remlpxem_new(y=yw, X=Xw, Z=Zw, dd.start=4, sigma2_e.start=8, maxiter=50, pxem=TRUE)

#-------------------------------------------
# average information algorithm using asreml
#-------------------------------------------
# asreml is commercial software available through VSN International (www.vsni.co.uk)
# asreml version 3 used here
# parameter updates (and convergence checking) performed manually outside of asreml so algorithm is the same as that presented by Gilmour et. al. 1995, i.e., without a rescaling of the average information matrix
library(asreml)

# Lamb weights
maxiter <- 50
lamb.asr.sv.sigma <- asreml(weight~1+DAMage+Line, random=~Sire, rcov=~idv(units), family=asreml.gaussian(dispersion=1), start.values=TRUE, data=lamb.df)
gam <- lamb.asr.sv.sigma$gammas.table
# SV1
gam$Value <- c(2,1,2)
# SV2
#gam$Value <- c(3,1,2)
AIupdate <- Kappa <- list()
Kappa.change <- NULL
Kappa[[1]] <- c(gam$Value[1], gam$Value[3])
for(w in 1:maxiter)
{
  lamb.asr.sigma <- asreml(weight~1+DAMage+Line, random=~Sire, rcov=~idv(units), family=asreml.gaussian(dispersion=1), G.param=gam, R.param=gam, control=asreml.control(maxiter=1), data=lamb.df)
  myAI <- diag(3)
  myAI[upper.tri(myAI, diag=TRUE)] <- lamb.asr.sigma$ai
  myAI <- myAI + t(myAI) - diag(diag(myAI)) 
  AIupdate[[w]] <- gam$Value + myAI%*%lamb.asr.sigma$score
  gam$Value <- AIupdate[[w]]
  Kappa[[w+1]] <- c(gam$Value[1], gam$Value[3])
  Kappa.change[w] <- sqrt( (t(Kappa[[w+1]] - Kappa[[w]])%*%(Kappa[[w+1]] - Kappa[[w]])) / (t(Kappa[[w]])%*%Kappa[[w]]) )
  if(Kappa.change[w] < 10^-8)
  {
    cat("Converged in",w,"iterations")
    break
  }
}

# Soybean variety trial yield example
maxiter <- 50
soy.asr.sv.sigma <- asreml(yield~1+gen, random=~block, rcov=~idv(units), family=asreml.gaussian(dispersion=1), start.values=TRUE, data=weiss.incblock)
gam <- soy.asr.sv.sigma$gammas.table
# SV1
gam$Value <- c(1,1,1)
# SV2
#gam$Value <- c(4,1,8)
AIupdate <- Kappa <- list()
Kappa.change <- NULL
Kappa[[1]] <- c(gam$Value[1], gam$Value[3])
for(w in 1:maxiter)
{
  soy.asr.sigma <- asreml(yield~1+gen, random=~block, rcov=~idv(units), family=asreml.gaussian(dispersion=1), G.param=gam, R.param=gam, control=asreml.control(maxiter=1), data=weiss.incblock)
  myAI <- diag(3)
  myAI[upper.tri(myAI, diag=TRUE)] <- soy.asr.sigma$ai
  myAI <- myAI + t(myAI) - diag(diag(myAI)) 
  AIupdate[[w]] <- gam$Value + myAI%*%soy.asr.sigma$score
  gam$Value <- AIupdate[[w]]
  Kappa[[w+1]] <- c(gam$Value[1], gam$Value[3])
  Kappa.change[w] <- sqrt( (t(Kappa[[w+1]] - Kappa[[w]])%*%(Kappa[[w+1]] - Kappa[[w]])) / (t(Kappa[[w]])%*%Kappa[[w]]) )
  if(Kappa.change[w] < 10^-8)
  {
    cat("Converged in",w,"iterations")
    break
  }
}


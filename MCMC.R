rCMLG <- function(H=matrix(rnorm(6),3), alpha=c(1,1,1), kappa=c(1,1,1)){
  ## This function will simulate from the cMLG distribution
  notValid <- T
  while(notValid){
    m <- length(kappa)
    w <- log(rgamma(m, shape=alpha, rate=kappa))
    HtH <- t(H)%*%H
    if(isDiagonal(HtH)) HtH <- Diagonal(x=diag(HtH))
    invHtH <- solve(HtH)
    tHw <- t(H)%*%w
    cand <- as.numeric(invHtH%*%tHw)
    if(!any(is.na(cand))) notValid=F
  }
  return(cand)
}

# dlMLG <- function(q, V, a, k){
#   det <- as.numeric(determinant(V%*%t(V), logarithm=T))
#   det <- det[1]*det[2]
#   ll <- -0.5*(det)+sum(a*log(k)) - sum(lgamma(a)) + 
#     t(a)%*%solve(V)%*%(q)-t(k)%*%exp(solve(V)%*%q)
#   return(ll)
# }

dlMLG <- function(q, V, a, k){
  V <- Matrix(V)
  Vinv <- solve((V))
  det <- as.numeric(determinant(V%*%t(V), logarithm=T))
  det <- det[1]*det[2]
  ll <- -0.5*(det)+sum(a*log(k)) - sum(lgamma(a)) + 
    t(a)%*%Vinv%*%(q)-t(k)%*%exp(Vinv%*%q)
  return(as.numeric(ll))
}


rTruncCMLG <- function(H=matrix(rnorm(6),3), alpha=c(1,1,1), kappa=c(1,1,1), cut=0 ){
  ## This function simulates from a cMLG distribution truncated from below
  repeat {
    temp <- rCMLG(H=H, alpha=alpha, kappa=kappa)
    if(all(temp > cut)) break
  }
  return(temp)
}

model1 <- function(y, X, Psi, wgt, iter=1000){
  n <- length(y)
  p <- ncol(X)
  r <- ncol(Psi)
  WGT <- Diagonal(n, wgt)
  XWX <- t(cbind(X, Psi))%*%WGT%*%cbind(X,Psi)
  XWY <- t(cbind(X, Psi))%*%WGT%*%y
  sig2<- 1
  sig2RE <- 1
  sig2Out <- rep(NA, iter)
  betaOut <- matrix(NA, nrow=iter, ncol=p+r)
  #pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    A <- XWX/sig2 + bdiag(Diagonal(p, 1/1000), Diagonal(r, 1/sig2RE))
    Ainv <- solve(A)
    b <- XWY/sig2
    beta <- betaOut[i,] <- as.vector(rmvnorm(1,Ainv%*%b,as.matrix(Ainv)))
    
    sig2 <- sig2Out[i] <- 1/rgamma(1,
                                   0.1 + n/2,
                                   0.1 + 0.5*crossprod(wgt*(y-cbind(X, Psi)%*%beta), (y-cbind(X, Psi)%*%beta)))
    sig2RE <- 1/rgamma(1,
                          0.1 + r/2,
                          0.1 + 0.5*sum(beta[-c(1:p)]^2))
    #setTxtProgressBar(pb, i)
  }
  return(list(Beta=betaOut, sig2=sig2Out))
}

# model1 <- function(y, X, wgt, iter=1000){
#   n <- length(y)
#   p <- ncol(X)
#   WGT <- Diagonal(n, wgt)
#   XWX <- t(X)%*%WGT%*%X
#   XWY <- t(X)%*%WGT%*%y
#   sig2<- 1
#   sig2Out <- rep(NA, iter)
#   betaOut <- matrix(NA, nrow=iter, ncol=p)
#   #pb <- txtProgressBar(min=0, max=iter, style=3)
#   for(i in 1:iter){
#     A <- XWX/sig2 + Diagonal(p, 1/1000)
#     Ainv <- solve(A)
#     b <- XWY/sig2
#     beta <- betaOut[i,] <- as.vector(rmvnorm(1,Ainv%*%b,as.matrix(Ainv)))
#     
#     sig2 <- sig2Out[i] <- 1/rgamma(1,
#                                    0.1 + n/2,
#                                    0.1 + 0.5*crossprod(wgt*(y-X%*%beta), (y-X%*%beta)))
#     #setTxtProgressBar(pb, i)
#   }
#   return(list(Beta=betaOut, sig2=sig2Out))
# }

model2 <- function(y, X, wgt, iter=1000, alpha=1000){
  n <- length(y)
  p <- ncol(X)
  sig2<- rep(1, n)
  sig2Out <- rep(NA, iter)
  betaOut <- etaOut <- matrix(NA, nrow=iter, ncol=p)
  #pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    WGT <- Diagonal(n, wgt/sig2)
    XWX <- t(X)%*%WGT%*%X
    XWY <- t(X)%*%WGT%*%y
    
    A <- XWX + Diagonal(p, 1/1000)
    Ainv <- solve(A)
    b <- XWY
    beta <- betaOut[i,] <- as.vector(rmvnorm(1,Ainv%*%b,as.matrix(Ainv)))
    
    H <- rbind(X, alpha^(-0.5)/1000*Diagonal(p))
    a <- c(wgt/2, rep(alpha,p))
    k <- c(wgt*(y-X%*%beta)^2/2, rep(alpha,p))
    etaOut[i,] <- eta <- rCMLG(H, a, k)
    
    sig2 <- exp(-X%*%eta)
    
    #setTxtProgressBar(pb, i)
  }
  return(list(Beta=betaOut, Eta=etaOut))
}


bvm2 <- function(y, X1, X2, psi1, psi2,  wgt=NULL, iter=1000, alpha=1000, burn=500, sig2b=1, w=1000, rho=1000){
  n <- length(y)
  if(is.null(wgt)) wgt <- rep(1, n)
  p <- ncol(X1)
  p2 <- ncol(X2)
  r <- ncol(psi1)
  r2 <- ncol(psi2)
  aB <- c(rep(alpha, r2), w)
  kB <- c(rep(alpha, r2), rho)
  xi1 <- cbind(X1, psi1)
  xi2 <- cbind(X2, psi2)
  sig2<- rep(1, n)
  sig2Out <- lppd <- matrix(NA, nrow=iter, ncol=n)
  sig2e1 <- sig2e2 <- 1
  betaOut <- matrix(NA, nrow=iter, ncol=p+r)
  etaOut <- matrix(NA, nrow=iter, ncol=p2+r2)
  LL <- rep(NA, iter)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    WGT <- Diagonal(n, wgt/sig2)
    XWX <- t(xi1)%*%WGT%*%xi1
    XWY <- t(xi1)%*%WGT%*%y
    
    bPCov <- bdiag(Diagonal(p, sig2b), Diagonal(r, sig2e1))
    A <- XWX + solve(bPCov)
    Ainv <- solve(A)
    b <- XWY
    beta <- betaOut[i,] <- as.vector(rmvnorm(1,Ainv%*%b,as.matrix(Ainv)))
    
    H <- rbind(wgt*xi2, alpha^(-0.5)*bdiag(Diagonal(p2, 1/sig2b), Diagonal(r2, 1/sig2e2)))
    a <- c(wgt*0+1/2, rep(alpha,p2+r2))
    k <- c(as.numeric(wgt*(y-xi1%*%beta)^2/2), rep(alpha,p2+r2))
    etaOut[i,] <- eta <- rCMLG( H, a, k)
    
    sig2 <-  sig2Out[i, ] <- as.numeric(exp(-xi2%*%eta))
    sig2e1 <- 1/rgamma(1, 0.5 + r/2, 0.5 + sum(beta[-c(1:p)]^2)/2)
    H <- matrix(c(alpha^(-0.5)*eta[-c(1:p2)], alpha^(-0.5)),ncol=1)
    sig2e2 <- rTruncCMLG(H=H, alpha=aB, kappa=kB)
    
    LL[i] <- sum(dnorm(y, mean=as.numeric(xi1%*%beta), sd=sqrt(sig2), log=T))
    lppd[i,] <- dnorm(y, mean=as.numeric(xi1%*%beta), sd=sqrt(sig2), log=F)
    
    setTxtProgressBar(pb, i)
  }
  pdic <- 2*(sum(dnorm(y, mean=as.numeric(xi1%*%colMeans(betaOut[-c(1:burn), ])), sd=sqrt(colMeans(sig2Out[-c(1:burn), ])), log=T))
             - mean(LL[-c(1:burn)]))
  DIC <- -2*sum(dnorm(y, mean=as.numeric(xi1%*%colMeans(betaOut[-c(1:burn), ])), sd=sqrt(colMeans(sig2Out[-c(1:burn), ])), log=T)) + 2*pdic
  WAIC <- sum(log(colMeans(lppd))) - sum(apply(log(lppd), 2, var))
  return(list(Beta=betaOut[-c(1:burn),], Eta=etaOut[-c(1:burn),], sig2=sig2Out[-c(1:burn),], DIC=DIC, WAIC=WAIC))
}

lbvm <- function(y, X1, X2, wgt=NULL, iter=1000, alpha=1000, burn=500, sig2b=1){
  n <- length(y)
  if(is.null(wgt)) wgt <- rep(1, n)
  p <- ncol(X1)
  p2 <- ncol(X2)
  sig2<- rep(1, n)
  eta <- rep(1, p2)
  beta <- rep(1, p)
  sig2Out <- lppd <- matrix(NA, nrow=iter, ncol=n)
  betaOut <- matrix(NA, nrow=iter, ncol=p)
  etaOut <- matrix(NA, nrow=iter, ncol=p2)
  LL <- rep(NA, iter)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    sinv <- rinvgauss(n, mean=sqrt(2/(exp(-X2%*%eta)*wgt*(y-X1%*%beta)^2)), shape=wgt*2/exp(-X2%*%eta))
    
    
    WGT <- Diagonal(n, wgt*sinv)
    XWX <- t(X1)%*%WGT%*%X1
    XWY <- t(X1)%*%WGT%*%y
    
    A <- XWX + Diagonal(p, 1/1000)
    Ainv <- forceSymmetric(solve(A))
    b <- XWY
    beta <- betaOut[i,] <- as.vector(rmvnorm(1,Ainv%*%b,as.matrix(Ainv)))
    
    H <- rbind(X2, alpha^(-0.5)/sig2b*Diagonal(p2))
    a <- c(wgt*0+1, rep(alpha,p2))
    k <- c(as.numeric(1/sinv), rep(alpha,p2))
    etaOut[i,] <- eta <- rCMLG( H, a, k)
    
    sig2 <-  sig2Out[i, ] <- as.numeric(exp(-X2%*%eta))
    
    LL[i] <- sum(dlaplace(y, location=as.numeric(X1%*%beta), scale=sqrt(sig2/2), log=T))
    lppd[i,] <- dlaplace(y, location=as.numeric(X1%*%beta), scale=sqrt(sig2/2))
    
    setTxtProgressBar(pb, i)
  }
  pdic <- 2*(sum(dlaplace(y, location=as.numeric(X1%*%colMeans(betaOut[-c(1:burn), ])), scale=colMeans(sqrt(sig2Out[-c(1:burn), ]/2)), log=T))
             - mean(LL[-c(1:burn)]))
  DIC <- -2*sum(dlaplace(y, location=as.numeric(X1%*%colMeans(betaOut[-c(1:burn), ])), scale=colMeans(sqrt(sig2Out[-c(1:burn), ]/2)), log=T)) + 2*pdic
  WAIC <- sum(log(colMeans(lppd))) - sum(apply(log(lppd), 2, var))
  return(list(Beta=betaOut[-c(1:burn),], Eta=etaOut[-c(1:burn),], sig2=sig2Out[-c(1:burn),], DIC=DIC, WAIC=WAIC))
}


FH_Fit <- function(Y, X, S2, sig2b=1000, iter=1000, burn=500){
  n <- length(Y)
  p <- ncol(X)
  tau2 <- 1
  tau2Out <- rep(NA, iter)
  beta1 <- rnorm(p, sd=100)
  beta1Out <- matrix(NA, nrow=p, ncol=iter)
  XtX <- t(X)%*%X
  thetaOut <- matrix(NA, nrow=n, ncol=iter)
  ll <- rep(NA, iter)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    thetaVar <- 1/(1/S2 + 1/tau2)
    thetaMu <- thetaVar*(Y/S2 + X%*%beta1/tau2)
    theta <- thetaOut[,i] <- rnorm(n, mean=thetaMu, sd=sqrt(thetaVar))
    
    betaVar <- solve(XtX/tau2 + Diagonal(p)/sig2b)
    betaMu <- as.numeric(betaVar%*%t(X)%*%theta/tau2)
    beta1 <- beta1Out[,i] <- as.numeric(rmvnorm(1, mean=betaMu, sigma=as.matrix(betaVar)))
    
    tau2 <- tau2Out[i] <- 1/rgamma(1,
                                   0.1 + n/2,
                                   0.1 + t(theta-X%*%beta1)%*%(theta-X%*%beta1)/2)
    ll[i] <- -2*sum(dnorm(Y, mean = theta, sd=sqrt(S2), log=T))
    
    setTxtProgressBar(pb, i)
  }
  DIC <- 2*mean(ll[-c(1:burn)]) + 2*sum(dnorm(Y, mean=rowMeans(thetaOut[,-c(1:burn)]), sd=sqrt(S2), log=T))
  return(list(Preds=thetaOut[,-c(1:burn)], Tau2=tau2Out[-c(1:burn)], Beta=beta1Out[,-c(1:burn)], DIC=DIC))
}

FH_Fit2 <- function(Y, X, S2, sig2b=1000, iter=1000, burn=500){
  n <- length(Y)
  p <- ncol(X)
  r <- n
  tau2 <- 1
  tau2Out <- rep(NA, iter)
  beta1Out <- matrix(NA, nrow=p, ncol=iter)
  eta1Out <- matrix(NA, nrow=r, ncol=iter)
  XP <- cbind(X, diag(n))
  thetaOut <- matrix(NA, nrow=n, ncol=iter)
  ll <- rep(NA, iter)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    ## sample mean function coefs
    var <- solve(t(XP)%*%Diagonal(x=1/S2)%*%XP + bdiag(1/sig2b*Diagonal(p), 1/tau2*Diagonal(r)))
    mu <- var%*%t(XP)%*%Diagonal(x=1/S2)%*%Y
    tmp <- as.numeric(rmvnorm(1, mean=mu, sigma=as.matrix(var)))
    beta1 <- beta1Out[,i] <- tmp[1:p]
    eta1 <- eta1Out[,i] <- tmp[-c(1:p)]
    theta <- thetaOut[,i] <- XP%*%c(beta1,eta1)
    
    tau2 <- tau2Out[i] <- 1/rgamma(1,
                                   0.1 + n/2,
                                   0.1 + t(theta-X%*%beta1)%*%(theta-X%*%beta1)/2)
    ll[i] <- -2*sum(dnorm(Y, mean = theta, sd=sqrt(S2), log=T))
    
    setTxtProgressBar(pb, i)
  }
  DIC <- 2*mean(ll[-c(1:burn)]) + 2*sum(dnorm(Y, mean=rowMeans(thetaOut[,-c(1:burn)]), sd=sqrt(S2), log=T))
  return(list(Preds=thetaOut[,-c(1:burn)], Tau2=tau2Out[-c(1:burn)], Beta=beta1Out[,-c(1:burn)], DIC=DIC))
}

saeFit <- function(Y, X, Psi, S2, SS, a=0.5, b=0.5, w=1000, rho=1000,sig2b=1000, iter=1000, burn=500, alpha=1000, propSD=1){
  df <- (SS-1)
  n <- length(Y)
  p <- ncol(X)
  r <- ncol(Psi)
  sig2 <- S2
  sig2e1 <- sig2e2 <- 0.2
  XP <- cbind(X, Psi)
  ICAR <- Diagonal(n)
  cholICAR <- solve(t(chol(ICAR)))
  aE <- c(rep(alpha, n), w)
  kE <- c(rep(alpha, n), rho)
  
  sig2Out <- matrix(NA, nrow=n, ncol=iter)
  thetaOut <- matrix(NA, nrow=n, ncol=iter)
  beta1Out <- beta2Out <- matrix(NA, nrow=p, ncol=iter)
  eta1Out <- eta2Out <- matrix(NA, nrow=r, ncol=iter)
  sig2e2Out <- rep(NA, iter)
  ll <- rep(NA, iter)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    ## sample mean function coefs
    var <- solve(t(XP)%*%Diagonal(x=1/sig2)%*%XP + bdiag(1/sig2b*Diagonal(p), 1/sig2e1*Diagonal(r)))
    mu <- var%*%t(XP)%*%Diagonal(x=1/sig2)%*%Y
    tmp <- as.numeric(rmvnorm(1, mean=mu, sigma=as.matrix(var)))
    beta1 <- beta1Out[,i] <- tmp[1:p]
    eta1 <- eta1Out[,i] <- tmp[-c(1:p)]
    theta <- thetaOut[,i] <- XP%*%c(beta1,eta1)
    
    
    ## sample var function coefs
    H <- rbind(XP, XP, bdiag(alpha^(-0.5)/sqrt(sig2b)*Diagonal(p),  alpha^(-0.5)/sqrt(sig2e2)*Diagonal(r)))
    a <- c(rep(1/2,n), df/2, rep(alpha,p+r))
    k <- c((Y-theta)^2/2, 1/2 *(SS-1)*S2, rep(alpha,p+r))
    tmp <- rCMLG(H, a, k)
    beta2 <- beta2Out[,i] <- tmp[1:p]
    eta2 <- eta2Out[,i] <- tmp[-c(1:p)]
    
    sig2 <- sig2Out[,i] <- exp(-XP%*%c(beta2, eta2))
    
    ## sample RE variances
    sig2e1 <- 1/rgamma(1, a+r/2, b+sum(eta1^2)/2)
    # H <- matrix(c(alpha^(-0.5)*eta2, 1),ncol=1)
    # sig2e2 <- sig2e2Out[i] <-  (1/rTruncCMLG(H=H, alpha=aE, kappa=kE))^2
    prop <- rtruncnorm(1, 0, mean=sig2e2, sd=propSD)
    ratio <- dlMLG(q=eta2, V=sqrt(alpha)*sqrt(prop)*diag(r), a=rep(alpha,r), k=rep(alpha,r)) -
      dlMLG(eta2, sqrt(alpha)*sqrt(sig2e2)*diag(r), rep(alpha,r), rep(alpha,r)) +
      dtruncnorm(sqrt(prop), 0, sd=5) - dtruncnorm(sqrt(sig2e2), 0, sd=5) +
      dtruncnorm(sig2e2, 0, mean=prop, sd=propSD) - dtruncnorm(prop, 0, mean=sig2e2, sd=propSD)
    if(runif(1) < exp(ratio)) sig2e2 <- prop
    sig2e2Out[i] <- sig2e2
    
    
    
    
    
    ll[i] <- -2*sum(dnorm(Y, mean = theta, sd=sqrt(sig2), log=T))
    setTxtProgressBar(pb, i)
  }
  DIC <- 2*mean(ll[-c(1:burn)]) + 2*sum(dnorm(Y, mean=rowMeans(thetaOut[,-c(1:burn)]), sd=sqrt(rowMeans(sig2Out[,-c(1:burn)])), log=T))
  return(list(Preds=thetaOut, DIC=DIC, sig2=sig2Out, beta1=beta1Out, beta2=beta2Out, eta2=eta2Out, sig2e2=sig2e2Out))
}



saeICAR <- function(Y, X, Psi, S2, SS, W, a=0.5, b=0.5, w=1, rho=1, sig2b=1000, iter=1000, burn=500, alpha=1000, propSD=1){
  n <- length(Y)
  p <- ncol(X)
  r <- ncol(Psi)
  sig2 <- S2
  sig2e1 <- 0.2
  sig2e2 <- 0.2
  XP <- cbind(X, Psi)
  dep <- 0.99
  ICAR <- diag(rowSums(W)) - dep*W
  cholICAR <- t(chol(ICAR))
  Vinv <- solve(cholICAR)
  aE <- c(rep(alpha, n), w)
  kE <- c(rep(alpha, n), rho)
  
  depOut <- rep(NA, iter)
  beta1 <- rep(1,p)
  sig2Out <- matrix(NA, nrow=n, ncol=iter)
  thetaOut <- matrix(NA, nrow=n, ncol=iter)
  beta1Out <- beta2Out <- matrix(NA, nrow=p, ncol=iter)
  eta1Out <- eta2Out <- matrix(NA, nrow=r, ncol=iter)
  sig2e2Out <- rep(NA,iter)
  ll <- rep(NA, iter)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    ## sample mean function coefs
    var <- solve(t(XP)%*%Diagonal(x=1/sig2)%*%XP + bdiag(1/sig2b*Diagonal(p), 1/sig2e1*ICAR))
    mu <- var%*%t(XP)%*%Diagonal(x=1/sig2)%*%Y
    tmp <- as.numeric(rmvnorm(1, mean=mu, sigma=as.matrix(var)))
    beta1 <- beta1Out[,i] <- tmp[1:p]
    eta1 <- eta1Out[,i] <- tmp[-c(1:p)]
    theta <- thetaOut[,i] <- XP%*%c(beta1,eta1)
    
    
    
    prop <- runif(1)
    ratio <- dmvnorm(eta1, sigma=sig2e1*solve(diag(rowSums(W)) - prop*W), log=T) -
      dmvnorm(eta1, sigma=sig2e1*solve(diag(rowSums(W)) - dep*W), log=T)
    if(runif(1) < exp(ratio)) dep <- prop
    depOut[i] <- dep
    ICAR <- diag(rowSums(W)) - dep*W
    
    
    ## sample var function coefs
    H <- rbind(XP, XP, bdiag(alpha^(-0.5)/sqrt(sig2b)*Diagonal(p),  alpha^(-0.5)/sqrt(sig2e2)*Diagonal(r)))
    a <- c(rep(1/2,n), (SS -1)/2, rep(alpha,p+r))
    k <- c((Y-theta)^2/2, 1/2 *(SS-1)*S2, rep(alpha,p+r))
    tmp <- rCMLG(H, a, k)
    beta2 <- beta2Out[,i] <- tmp[1:p]
    eta2 <- eta2Out[,i] <- tmp[-c(1:p)]
    
    sig2 <- sig2Out[,i] <- exp(-XP%*%c(beta2, eta2))
    
    ## sample RE variances
    sig2e1 <- 1/rgamma(1, a+r/2, b+t(eta1)%*%ICAR%*%eta1/2)
    # H <- matrix(c(alpha^(-0.5)*eta2, 1),ncol=1)
    # sig2e2 <- sig2e2Out[i] <-  (1/rTruncCMLG(H=H, alpha=aE, kappa=kE))^2
    prop <- rtruncnorm(1, 0, mean=sig2e2, sd=propSD)
    ratio <- dlMLG(eta2, sqrt(alpha)*sqrt(prop)*diag(r), rep(alpha,r), rep(alpha,r)) -
      dlMLG(eta2, sqrt(alpha)*sqrt(sig2e2)*diag(r), rep(alpha,r), rep(alpha,r)) +
      dtruncnorm(sqrt(prop), 0, sd=5) - dtruncnorm(sqrt(sig2e2), 0, sd=5) +
      dtruncnorm(sig2e2, 0, mean=prop, sd=propSD) - dtruncnorm(prop, 0, mean=sig2e2, sd=propSD)
    if(runif(1) < exp(ratio)) sig2e2 <- prop
    sig2e2Out[i] <- sig2e2
    
    
    
    
    
    ll[i] <- -2*sum(dnorm(Y, mean = theta, sd=sqrt(sig2), log=T))
    setTxtProgressBar(pb, i)
  }
  DIC <- 2*mean(ll[-c(1:burn)]) + 2*sum(dnorm(Y, mean=rowMeans(thetaOut[,-c(1:burn)]), sd=sqrt(rowMeans(sig2Out[,-c(1:burn)])), log=T))
  return(list(Preds=thetaOut, DIC=DIC, sig2=sig2Out, beta1=beta1Out, beta2=beta2Out, eta1=eta1Out, eta2=eta2Out, dep=depOut))
}


saeICAR2 <- function(Y, X, Psi, S2, SS, W, a=0.5, b=0.5, w=1, rho=1, sig2b=1000, iter=1000, burn=500, alpha=1000, propSD=1){
  df <- (SS-1)
  n <- length(Y)
  p <- ncol(X)
  r <- ncol(Psi)
  sig2 <- S2
  sig2e1 <- 0.2
  sig2e2 <- 0.2
  XP <- cbind(X, Psi)
  #ICAR <- chol(diag(rowSums(W)) - .99*W)
  ICAR <- (diag(rowSums(W)) - W)
  cholICAR <- abs(sqrtm(ICAR))
  Vinv <- solve(cholICAR)
  aE <- c(rep(alpha, n), w)
  kE <- c(rep(alpha, n), rho)
  
  depOut <- rep(NA, iter)
  beta1 <- rep(1,p)
  sig2Out <- matrix(NA, nrow=n, ncol=iter)
  thetaOut <- matrix(NA, nrow=n, ncol=iter)
  beta1Out <- beta2Out <- matrix(NA, nrow=p, ncol=iter)
  eta1Out <- eta2Out <- matrix(NA, nrow=r, ncol=iter)
  sig2e2Out <- rep(NA,iter)
  ll <- rep(NA, iter)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    ## sample mean function coefs
    # t1 <- Sys.time()
    # var <- as.matrix(forceSymmetric(solve(t(XP)%*%Diagonal(x=1/sig2)%*%XP + bdiag(1/sig2b*Diagonal(p), 1/sig2e1*ICAR))))
    # mu <- var%*%t(XP)%*%Diagonal(x=1/sig2)%*%Y
    # tmp <- as.numeric(rmvnorm(1, mean=mu, sigma=as.matrix(var)))
    # beta1 <- beta1Out[,i] <- tmp[1:p]
    # eta1 <- eta1Out[,i] <- tmp[-c(1:p)]
    # theta <- thetaOut[,i] <- XP%*%c(beta1,eta1)
    # Sys.time()-t1
    
    #t1 <- Sys.time()
    Q <- forceSymmetric(t(XP)%*%Diagonal(x=1/sig2)%*%XP + bdiag(1/sig2b*Diagonal(p), 1/sig2e1*ICAR))
    U <- chol(Q)
    a <- t(XP)%*%Diagonal(x=1/sig2)%*%Y
    bb <- rnorm(nrow(Q))
    tmp <- backsolve(U, backsolve(U, a, transpose = TRUE) + bb)
    beta1 <- beta1Out[,i] <- tmp[1:p]
    eta1 <- eta1Out[,i] <- tmp[-c(1:p)]
    theta <- thetaOut[,i] <- XP%*%c(beta1,eta1)
    #Sys.time()-t1
    
    
    
    
    
    ## sample var function coefs
    H <- rbind(XP, XP, bdiag(alpha^(-0.5)/sqrt(sig2b)*Diagonal(p),  alpha^(-0.5)/sqrt(sig2e2)*t(cholICAR)))
    a <- c(rep(1/2,n), df/2, rep(alpha,p+r))
    k <- c((Y-theta)^2/2, 1/2 *(SS-1)*S2, rep(alpha,p+r))
    tmp <- rCMLG(H, a, k)
    beta2 <- beta2Out[,i] <- tmp[1:p]
    eta2 <- eta2Out[,i] <- tmp[-c(1:p)]

    
    sig2 <- sig2Out[,i] <- exp(-XP%*%c(beta2, eta2))
    
    ## sample RE variances
    sig2e1 <- 1/rgamma(1, a+r/2, b+t(eta1)%*%ICAR%*%eta1/2)
    # H <- matrix(c(alpha^(-0.5)*eta2, 1),ncol=1)
    # sig2e2 <- sig2e2Out[i] <-  (1/rTruncCMLG(H=H, alpha=aE, kappa=kE))^2
    prop <- rtruncnorm(1, 0, mean=sig2e2, sd=propSD)
    ratio <- dlMLG(eta2, sqrt(alpha)*sqrt(prop)*Vinv, rep(alpha,r), rep(alpha,r)) -
      dlMLG(eta2, sqrt(alpha)*sqrt(sig2e2)*Vinv, rep(alpha,r), rep(alpha,r)) +
      dtruncnorm(sqrt(prop), 0, sd=5) - dtruncnorm(sqrt(sig2e2), 0, sd=5) +
      dtruncnorm(sig2e2, 0, mean=prop, sd=propSD) - dtruncnorm(prop, 0, mean=sig2e2, sd=propSD)
    if(runif(1) < exp(ratio)) sig2e2 <- prop
    sig2e2Out[i] <- sig2e2
    
    
    
    
    
    ll[i] <- -2*sum(dnorm(Y, mean = theta, sd=sqrt(sig2), log=T))
    setTxtProgressBar(pb, i)
  }
  DIC <- 2*mean(ll[-c(1:burn)]) + 2*sum(dnorm(Y, mean=rowMeans(thetaOut[,-c(1:burn)]), sd=sqrt(rowMeans(sig2Out[,-c(1:burn)])), log=T))
  return(list(Preds=thetaOut, DIC=DIC, sig2=sig2Out, beta1=beta1Out, beta2=beta2Out, eta1=eta1Out, eta2=eta2Out, dep=depOut, sig2e2=sig2e2Out))
}

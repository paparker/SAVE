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



saeICAR2 <- function(Y, X, Psi, S2, SS, W, a=0.5, b=0.5, w=1, rho=1, sig2b=1000, iter=1000, burn=500, alpha=1000, propSD=1){
  df <- (SS-1)
  n <- length(Y)
  p <- ncol(X)
  r <- ncol(Psi)
  sig2 <- S2
  sig2e1 <- 0.2
  sig2e2 <- 0.2
  XP <- cbind(X, Psi)
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

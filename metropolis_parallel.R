# metropolis fixed phi joint updating


library(MASS)
library(mvtnorm)
library(invgamma)
library(truncnorm)

ARGS <- commandArgs(trailingOnly = TRUE)
numSites <- as.numeric(ARGS[1])
numTotal <- as.numeric(ARGS[2])
a <- as.numeric(ARGS[3])
r <- as.numeric(ARGS[4])
sigmasq <- as.numeric(ARGS[5])
tausq <- as.numeric(ARGS[6])





### functions 

makeB <- function(B){  
  b1 <- B[1,1]
  b2 <- B[2,1]
  b3 <- B[1,2]
  b4 <- B[2,2]
  pwr <- Dx^2 * b1 + Dx * Dy * (b2 + b3) + Dy^2 * b4
  return(pwr)
}



# simulate data
par(mfrow=c(1,1))
set.seed(1127)

sites <- cbind(runif(numTotal, 0, 1), runif(numTotal, 0, 1))
hMat <- NULL
mu <- 0
phi <- 3/0.3
amin <- 1
amax <- amin*r
rotationMat <- matrix(c(cos(a),-sin(a),sin(a),cos(a)),nrow=2,ncol=2)
aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
A <- rotationMat %*% aMat
B <- A %*% t(A)
Sigma <- matrix(0,nrow=numTotal, ncol=numTotal)
for(m in 1:numTotal){
  for(n in 1:numTotal){
    h <- sites[m,] - sites[n,]
    Sigma[m,n] <- sigmasq * exp(-phi * sqrt(t(h) %*% B %*% h))
  }
}
Sigma <- Sigma + diag(tausq, nrow=numTotal, ncol=numTotal)
y <- mvrnorm(1, rep(0, numTotal), Sigma)
ind_test <- sample(1:numTotal, size=numTotal-numSites)
ind_train <- setdiff(1:numTotal, ind_test)
sites_test <- sites[ind_test,]
sites_train <- sites[ind_train,]
y_test <- y[ind_test]
y_train <- y[ind_train]

## distance matrix for x and y coordinate
Dx <- matrix(0, nrow=numSites, ncol=numSites)
Dy <- matrix(0, nrow=numSites, ncol=numSites)
for (m in 1:(numSites-1)) {
  for (n in (m+1):numSites) {
    h <- sites_train[m,] - sites_train[n,]
    Dx[m,n] <- h[1]
    Dx[n,m] <- h[1]
    Dy[m,n] <- h[2]
    Dy[n,m] <- h[2]
  }
}





# propose sigmasq

# initialize
numSim <- 100000
mcmc.sigma <- matrix(1, nrow = numSim, ncol = 1)
mcmc.a <- matrix(pi/3, nrow = numSim, ncol = 1)
mcmc.r <- matrix(4, nrow = numSim, ncol = 1)
mcmc.tau <- matrix(0.2, nrow = numSim, ncol = 1)
mcmc.mu <- matrix(0, nrow=numSim, ncol=numSites)
mcmc.phi <- matrix(3/0.3, nrow=numSim, ncol=1)



# Metropolis

for(i in 2:numSim){
  
  # update phi 
  var <- 2
  phi_star <- rlnorm(1, log(mcmc.phi[i-1,]), sqrt(var))
  amin <- 1
  amax <- amin*mcmc.r[i-1,]
  rotationMat <- matrix(c(cos(mcmc.a[i-1,]),-sin(mcmc.a[i-1,]),sin(mcmc.a[i-1,]),cos(mcmc.a[i-1,])),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  # covariances
  Sigma_star <- mcmc.sigma[i-1,] * exp(-phi_star * sqrt(makeB(B)))
  Sigma_prev <- mcmc.sigma[i-1,] * exp(-mcmc.phi[i-1,] * sqrt(makeB(B)))
  diag(Sigma_star) <- mcmc.sigma[i-1,] + mcmc.tau[i-1,]
  diag(Sigma_prev) <- mcmc.sigma[i-1,] + mcmc.tau[i-1,]
  # calulate ratio
  mvn_star <- dmvnorm(y_train, mcmc.mu[i-1,], Sigma_star, log=FALSE)
  mvn_prev <- dmvnorm(y_train, mcmc.mu[i-1,], Sigma_prev, log=FALSE)
  mvn_star <- dmvnorm(y_train, mcmc.mu[i-1,], Sigma_star)
  mvn_prev <- dmvnorm(y_train, mcmc.mu[i-1,], Sigma_prev)
  prior_star <- dunif(phi_star, 3/0.7, 3/0.2)
  prior_prev <- dunif(mcmc.phi[i-1,], 3/0.7, 3/0.2)
  ratio <- mvn_star * prior_star / (mvn_prev * prior_prev)
  u <- runif(1)
  if(log(u) < log(ratio)){
    mcmc.phi[i,] <- phi_star
  } else {
    mcmc.phi[i,] <- mcmc.phi[i-1,]
  }
  
  
  
  # update sigma sqaured
  
  var <- 3
  sigma_star <- rlnorm(1, log(mcmc.sigma[i-1,]), sqrt(var))
  #print("sigma")
  #print(c(sigma_star,mcmc.sigma[i-1,]))
  amin <- 1
  amax <- amin*mcmc.r[i-1,]
  rotationMat <- matrix(c(cos(mcmc.a[i-1,]),-sin(mcmc.a[i-1,]),sin(mcmc.a[i-1,]),cos(mcmc.a[i-1,])),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  # covariances
  Sigma_star <- sigma_star * exp(-mcmc.phi[i,] * sqrt(makeB(B)))
  Sigma_prev <- mcmc.sigma[i-1,] * exp(-mcmc.phi[i,] * sqrt(makeB(B)))
  diag(Sigma_star) <- sigma_star + mcmc.tau[i-1,]
  diag(Sigma_prev) <- mcmc.sigma[i-1,] + mcmc.tau[i-1,]
  
  # tan(theta) falls between -inf, to inf, as theta goes from 0 to pi, theta = arctan(z)
  # calulate ratio, log normal
  mvn_star <- dmvnorm(y_train, mcmc.mu[i-1,], Sigma_star)
  mvn_prev <- dmvnorm(y_train, mcmc.mu[i-1,], Sigma_prev)
  prior_star <- dinvgamma(sigma_star, shape=1, scale=1)
  prior_prev <- dinvgamma(mcmc.sigma[i-1,], shape=1, scale=1)
  ratio <- mvn_star * prior_star / (mvn_prev * prior_prev)
  u <- runif(1)
  if(log(u) < log(ratio)){
    mcmc.sigma[i,] <- sigma_star
  } else {
    mcmc.sigma[i,] <- mcmc.sigma[i-1,]
  }
  
  
  
  # update tau sqaured 
  var <- 3
  tau_star <- rlnorm(1, log(mcmc.tau[i-1,]), sqrt(var))
  #print("tau")
  #print(c(tau_star,mcmc.tau[i-1,]))
  amin <- 1
  amax <- amin*mcmc.r[i-1,]
  rotationMat <- matrix(c(cos(mcmc.a[i-1,]),-sin(mcmc.a[i-1,]),sin(mcmc.a[i-1,]),cos(mcmc.a[i-1,])),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  # covariance sampled
  Sigma <- mcmc.sigma[i,] * exp(-mcmc.phi[i,] * sqrt(makeB(B)))
  Sigma_prev <- Sigma
  Sigma_star <- Sigma
  diag(Sigma_star) <- mcmc.sigma[i,] + tau_star
  diag(Sigma_prev) <- mcmc.sigma[i,] + mcmc.tau[i-1,]
  # calulate ratio
  mvn_star <- dmvnorm(y_train, mcmc.mu[i-1,], Sigma_star, log=FALSE)
  mvn_prev <- dmvnorm(y_train, mcmc.mu[i-1,], Sigma_prev, log=FALSE)
  prior_star <- dinvgamma(tau_star, shape=1, scale=1)
  prior_prev <- dinvgamma(mcmc.tau[i-1,], shape=1, scale=1)
  ratio <- mvn_star * prior_star / (mvn_prev * prior_prev)
  u <- runif(1)
  if(log(u) < log(ratio)){
    mcmc.tau[i,] <- tau_star
  } else {
    mcmc.tau[i,] <- mcmc.tau[i-1,]
  }
  
  
  
  
  # update angle and ratio
  var1 <- 10
  tt <- rnorm(1, mcmc.a[i-1,], sqrt(var1))
  a_star <- tt%%(pi)
  var2 <- 8
  #r_star <- rlnorm(1,log(mcmc.r[i-1,]), sqrt(var2))
  r_star <- rtruncnorm(1, a=0, b=Inf, mean = mcmc.r[i-1,], sd = sqrt(var2))
  #print("a")
  #print(c(a_star,mcmc.a[i-1,]))
  #print("r")
  #print(c(r_star,mcmc.r[i-1,]))
  amin <- 1
  # angle ratio sampled
  amax <- amin*r_star
  rotationMat <- matrix(c(cos(a_star),-sin(a_star),sin(a_star),cos(a_star)),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B_star <- A %*% t(A)
  # angle previous
  amax <- amin*mcmc.r[i-1,]
  rotationMat <- matrix(c(cos(mcmc.a[i-1,]),-sin(mcmc.a[i-1,]),sin(mcmc.a[i-1,]),cos(mcmc.a[i-1,])),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B_prev <- A %*% t(A)
  #covariances
  Sigma_star <- mcmc.sigma[i,] * exp(-mcmc.phi[i,] * sqrt(makeB(B_star)))
  Sigma_prev <- mcmc.sigma[i,] * exp(-mcmc.phi[i,] * sqrt(makeB(B_prev)))
  diag(Sigma_star) <- mcmc.sigma[i,] + mcmc.tau[i,]
  diag(Sigma_prev) <- mcmc.sigma[i,] + mcmc.tau[i,]
  # calulate ratio
  mvn_star <- dmvnorm(y_train, mcmc.mu[i-1,], Sigma_star, log=TRUE)
  mvn_prev <- dmvnorm(y_train, mcmc.mu[i-1,], Sigma_prev, log=TRUE)
  prior_star <- dunif(a_star, 0, pi) * dinvgamma(r_star, shape=1, scale=1)
  prior_prev <- dunif(mcmc.a[i-1,], 0, pi) * dinvgamma(mcmc.r[i-1,], shape=1, scale=1)
  ratio <- exp(mvn_star-mvn_prev) * prior_star / prior_prev
  u <- runif(1)
  if(log(u) < log(ratio)){
    mcmc.a[i,] <- a_star
    mcmc.r[i,] <- r_star
  } else {
    mcmc.a[i,] <- mcmc.a[i-1,]
    mcmc.r[i,] <- mcmc.r[i-1,]
  }
  
  
  
  
  # update mu
  amax <- amin*mcmc.r[i,]
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  rotationMat <- matrix(c(cos(mcmc.a[i,]),-sin(mcmc.a[i,]),sin(mcmc.a[i,]),cos(mcmc.a[i,])),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B_prev <- A %*% t(A)
  #covariances
  Sigma_prev <- mcmc.sigma[i,] * exp(-mcmc.phi[i,] * sqrt(makeB(B_prev)))
  diag(Sigma_prev) <- mcmc.sigma[i,] + mcmc.tau[i,]
  var <- chol2inv(chol(diag(0.001, numSites, numSites) +  numSites * chol2inv(chol(Sigma_prev))))
  mean <- var %*% (chol2inv(chol(Sigma_prev)) %*% y_train)
  newmu <- mvrnorm(1, mean, var)
  mcmc.mu[i,] <- newmu
  
}

numBurn <- 80000
numKeep <- seq(numBurn+20, numSim, by=20)
a_samp <- mcmc.a[numKeep]
r_samp <- mcmc.r[numKeep]
sig_samp <- mcmc.sigmasq[numKeep]
tau_samp <- mcmc.tausq[numKeep]
phi_samp <- mcmc.phi[numKeep]
mu_samp <- mcmc.mu[numKeep,]



# marginalize out w
makePredictionAniso_v2 <- function(numSample, newSites){
  
  sigmasq <- sig_samp[numSample]
  tausq <- tau_samp[numSample]
  mu <- mu_samp[numSample,]
  a <- a_samp[numSample]
  r <- r_samp[numSample]
  phi <- phi_samp[numSample]
  amin <- 1
  amax <- amin*r
  rotationMat <- matrix(c(cos(a),-sin(a),sin(a),cos(a)),nrow=2,ncol=2)
  aMat <- matrix(c(1/amax,0,0,1/amin),nrow=2,ncol=2)
  A <- rotationMat %*% aMat
  B <- A %*% t(A)
  SigmaObs <- sigmasq * exp(-phi * sqrt(makeB(B)))
  diag(SigmaObs) <- sigmasq + tausq

  S <- solve(SigmaObs)
  
  predictions <- NULL
  for(j in 1:nrow(newSites)){
    newSite <- newSites[j,]
    distPO <- matrix(0, nrow=numSites, ncol = 2)
    for(i in 1:numSites){
      distPO[i,] <- newSite - sites_train[i,]
    }
    SigmaPO <- matrix(0, nrow=1, ncol=numSites)
    for(i in 1:numSites){
      SigmaPO[,i] <- sigmasq * exp(-phi * sqrt(distPO[i,] %*% B %*% distPO[i,]))
    }
    mean_yPred <- mu + SigmaPO %*% S %*% (y_train - rep(mu, numSites))
    var_yPred <- sigmasq + tausq - SigmaPO %*% S %*% t(SigmaPO)
    yPred <- rnorm(1, mean_yPred, sqrt(var_yPred))
    predictions <- c(predictions, yPred)
  }
  
  return(predictions)
}


# isotropy krigging 
makePredictionIso <- function(numSample, newSites){
  
  sigmasq <- sub.samps.iso[numSample,"sigmasq"]
  tausq <- sub.samps.iso[numSample,"tausq"]
  mu <- sub.samps.iso[numSample, "mu"]
  phi <- sub.samps.iso[numSample, "phi"]
  
  SigmaObs <- matrix(0,nrow=numSites, ncol=numSites)
  for (i in 1:(numSites-1)) {
    for (j in (i+1):numSites) {
      SigmaObs[i,j] <- sigmasq * exp(-phi * distX[i,j])
    }
  }
  
  for (i in 1:(numSites-1)) {
    for (j in (i+1):numSites) {
      SigmaObs[j,i] <- sigmasq * exp(-phi * distX[i,j])
    }
  }
  
  for(k in 1:numSites){
    SigmaObs[k,k] <- sigmasq + tausq
  }
  
  S <- solve(SigmaObs)
  
  predictions <- NULL
  for(j in 1:nrow(newSites)){
    newSite <- newSites[j,]
    distPO <- matrix(0, nrow=numSites, ncol = 1)
    for(i in 1:numSites){
      distPO[i,] <- dist(rbind(newSite,sites_train[i,]))
    }
    SigmaPO <- matrix(0, nrow=1, ncol=numSites)
    for(i in 1:numSites){
      SigmaPO[,i] <- sigmasq * exp(-phi * distPO[i,])
    }
    mean_yPred <- mu + SigmaPO %*% S %*% (y_train - rep(mu, numSites))
    var_yPred <- sigmasq + tausq - SigmaPO %*% S %*% t(SigmaPO)
    yPred <- rnorm(1, mean_yPred, sqrt(var_yPred))
    predictions <- c(predictions, yPred)
  }
  
  return(predictions)
  
}



nMCMC <- nrow(sub.samps.ani)
distX <- as.matrix(dist(sites_train))
predsIso <- t(sapply(1:nMCMC, makePredictionIso, newSites=sites_test))
predsAniso_v1 <- t(sapply(1:nMCMC, makePredictionAniso_v1, newSites=sites_test))
predsAniso_v2 <- t(sapply(1:nMCMC, makePredictionAniso_v2, newSites=sites_test))


save(predsIso, file="predsIso_ar_100.Rdata")
save(predsAniso_v1, file="predsAniso_ar_100_v1.Rdata")
save(predsAniso_v2, file="predsAniso_ar_100_v2.Rdata")


# Isotropy

# 1. empirical coverage
yObs <- y_test
numPred <- length(y_test)
qt <- apply(predsIso, 2, function(x) quantile(x, probs = c(0.05, 0.95)))
ave <- apply(predsIso, 2, mean)
empCov <- data.frame(cbind(ave, t(qt), yObs))
colnames(empCov) <- c("Mean", "Lower", "Higher", "True")
empCov$capture <- empCov$Lower <= empCov$True & empCov$Higher >= empCov$True 
pdf("empCovIso.pdf")
ggplot(empCov, aes(y=True, x=1:numPred, color=capture)) +
  xlab("index") +
  geom_errorbar(aes(ymax=Higher, ymin=Lower), width=0, color='black', alpha=0.3, size=2) +
  geom_point(size=5) +
  labs(x = "Index of New Sites", y="Predcited/Observed Value") +
  ggtitle(paste("Empirical Coverage of Isotropic Model =",sum(empCov$capture)/numPred))
dev.off()
ecIso <- sum(empCov$capture)/numPred


# 2. PMSE
mseIso <- mean((empCov$Mean-empCov$True)^2)

# 3. CRPS
crpsIso <- crps_test(empCov$Mean, empCov$True)

# Anisotropy

# 1. empirical coverage
qt <- apply(predsAniso_v2, 2, function(x) quantile(x, probs = c(0.05, 0.95)))
ave <- apply(predsAniso_v2, 2, mean)
empCov <- data.frame(cbind(ave, t(qt), yObs))
colnames(empCov) <- c("Mean", "Lower", "Higher", "True")
empCov$capture <- empCov$Lower <= empCov$True & empCov$Higher >= empCov$True 
pdf("empCovAniso.pdf")
ggplot(empCov, aes(y=True, x=1:length(y_test), color=capture)) +
  xlab("index") +
  geom_errorbar(aes(ymax=Higher, ymin=Lower), width=0, color='black', alpha=0.3, size=2) +
  geom_point(size=5) +
  labs(x = "Index of New Sites", y="Predcited/Observed Value") +
  ggtitle(paste("Empirical Coverage of Anisotropic Model =",sum(empCov$capture)/numPred))
dev.off()
ecAni <- sum(empCov$capture)/numPred


# 2. PMSE
mseAni <- mean((empCov$Mean-empCov$True)^2)

# 3. CRPS

crpsAni <- crps_test(empCov$Mean, empCov$True)

save(ecIso, mseIso, crpsIso, ecAni, mseAni, crpsAni, file="resPrediction_ar_1.Rdata")













name <- paste0(numSites,numSim,"aniso.Rdata")
save(mcmc.mu, mcmc.a, mcmc.r, mcmc.sigma, mcmc.tau, file=name)




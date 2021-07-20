#### This file contains all the functions related to the project ``zeroInf"

#### 1.  Global Parameters ####

ProjectName <- "zeroInf"

#### 2. Auxiliary Functions ####

cexpexp <- function(x){
  return(1-exp(-exp(x)))
}

expit <- function(x){
  value <- exp(x)/(1+exp(x))
  ifelse(is.na(value),1,value) 
}


#### 3. Main Functions ####

#### 3.1 Functions for Generating Simulation Data ####

### This function is for generating the simulated data of the FIRST simulation studies with binary covariates

datageneration <- function(seedid, nsample, beta1, beta2, alphaplus, alphaminus){
  
  set.seed(seed_i[seedid])
  
  ### Response model 
  
  X1 <- rbinom(nsample,1,0.5)
  X2 <- rbinom(nsample,1,0.5)
  Xplus <- rbinom(nsample,1,0.5)
  # Xminus <- rbinom(nsample,1,0.5)
  Xminus <- runif(nsample,0.5,1.5)
  
  covariates <- data.frame(intercept=rep(1,nsample),X1=X1, X2=X2, Xplus=Xplus, Xminus=Xminus)
  
  beta1 <- beta1
  beta2 <- beta2
  
  mu1 <- as.matrix(covariates[,c("intercept","X1")]) %*% t(t(beta1))
  mu2 <- as.matrix(covariates[,c("intercept","X2")]) %*% t(t(beta2))
  
  # UforY <- runif(nsample,0,1)
  
  U1 <- unlist(lapply(exp(mu1),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  U2 <- unlist(lapply(exp(mu2),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  Y <- (U1 > 0) *  U2
  
  ### Measurement error model
  
  # alphaplus <- alphaplus
  # alphaminus <- alphaminus
  
  muZplus <- as.matrix(covariates[,c("Xplus")]) * alphaplus 
  # muZminus <- as.matrix(covariates[,c("intercept","Xminus")]) %*% t(t(alphaminus))
  muZminus <- as.matrix(covariates[,c("Xminus")]) * alphaminus
  
  Zplus <- unlist(lapply(exp(muZplus),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  
  # muZminuspnorm <- pnorm(muZminus)
  muZminusexpit <- expit(muZminus)
  
  Zminus <- unlist(lapply(1:nsample,FUN=function(x){
    rbinom(1, size = Y[x], prob = muZminusexpit[x])
  }))
  
  Ystar <- Y + Zplus - Zminus
  
  return(data.frame(Ystar=Ystar,Y=Y,covariates,Zplus=Zplus,Zminus=Zminus,U1=U1,U2=U2))
}

### This function is for generating the simulated data of the FIRST simulation studies with continuous covariates

datagencontinue <- function(seedid, nsample, beta1, beta2, alphaplus, alphaminus){
  
  set.seed(seed_i[seedid])
  
  ### Response model 
  
  # X1 <- runif(nsample,0,1)
  X1 <- rbinom(nsample,1,0.5)
  X2 <- runif(nsample,0,1)
  Xplus <-  runif(nsample,0,1)
  # Xminus <- rbinom(nsample,1,0.5)
  Xminus <- runif(nsample,0,2)
  
  covariates <- data.frame(intercept=rep(1,nsample),X1=X1, X2=X2, Xplus=Xplus, Xminus=Xminus)
  
  beta1 <- beta1
  beta2 <- beta2
  
  mu1 <- as.matrix(covariates[,c("intercept","X1")]) %*% t(t(beta1))
  mu2 <- as.matrix(covariates[,c("intercept","X2")]) %*% t(t(beta2))
  
  # UforY <- runif(nsample,0,1)
  
  U1 <- unlist(lapply(exp(mu1),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  U2 <- unlist(lapply(exp(mu2),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  Y <- (U1 > 0) *  U2
  
  ### Measurement error model
  
  
  # muZplus <- as.matrix(covariates[,c("intercept","Xplus")]) %*% t(t(alphaplus))
  muZplus <- as.matrix(covariates[,c("Xplus")]) * alphaplus
  muZminus <- as.matrix(covariates[,c("Xminus")]) * alphaminus
  
  Zplus <- unlist(lapply(exp(muZplus),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  
  muZminuspnorm <- expit(muZminus)
  Zminus <- unlist(lapply(1:nsample,FUN=function(x){
    rbinom(1, size = Y[x], prob = muZminuspnorm[x])
  }))
  
  Ystar <- Y + Zplus - Zminus
  
  return(data.frame(Ystar=Ystar,Y=Y,covariates,Zplus=Zplus,Zminus=Zminus,U1=U1,U2=U2))
}

### This function is for generating the simulated data of the SECOND simulation studies with continuous covariates

dataSim3 <- function(seedid, nsample, beta1, beta2, alphaplus, alphaminus){
  
  set.seed(seed_i[seedid])
  
  ### Response model 
  
  # X1 <- runif(nsample,0,1)
  X1 <- rbinom(nsample,1,0.5)
  X2 <- runif(nsample,0,1)
  Xplus <-  runif(nsample,0,1)
  # Xminus <- rbinom(nsample,1,0.5)
  Xminus <- runif(nsample,0,2)
  
  covariates <- data.frame(intercept=rep(1,nsample),X1=X1, X2=X2, Xplus=Xplus, Xminus=Xminus)
  
  beta1 <- beta1
  beta2 <- beta2
  
  mu1 <- as.matrix(covariates[,c("intercept","X1")]) %*% t(t(beta1))
  mu2 <- as.matrix(covariates[,c("intercept","X2")]) %*% t(t(beta2))
  
  # UforY <- runif(nsample,0,1)
  
  U1 <- unlist(lapply(exp(mu1),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  U2 <- unlist(lapply(exp(mu2),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  Y <- (U1 > 0) *  U2
  
  ### Measurement error model
  
  
  muZplus <- as.matrix(covariates[,c("intercept","Xplus")]) %*% t(t(alphaplus))
  muZminus <- as.matrix(covariates[,c("intercept","Xminus")]) %*% t(t(alphaminus))
  # muZplus <- as.matrix(covariates[,c("Xplus")]) * alphaplus
  # muZminus <- as.matrix(covariates[,c("Xminus")]) * alphaminus
  
  Zplus <- unlist(lapply(exp(muZplus),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  
  muZminuspnorm <- expit(muZminus)
  Zminus <- unlist(lapply(1:nsample,FUN=function(x){
    rbinom(1, size = Y[x], prob = muZminuspnorm[x])
  }))
  
  Ystar <- Y + Zplus - Zminus
  
  return(data.frame(Ystar=Ystar,Y=Y,covariates,Zplus=Zplus,Zminus=Zminus,U1=U1,U2=U2))
}

### This function is for generating the simulated data with External Validation Sample

dataSim4 <- function(seedid, nmain, nvalid, beta1, beta2, alphaplus, alphaminus){
  
  nsample <- nmain + nvalid
  set.seed(seed_i[seedid])
  
  ### Response model 
  
  # X1 <- runif(nsample,0,1)
  X1 <- rbinom(nsample,1,0.5)
  X2 <- runif(nsample,0,1)
  
  covariates <- data.frame(intercept=rep(1,nsample), X1=X1, X2=X2)
  
  beta1 <- beta1
  beta2 <- beta2
  
  mu1 <- as.matrix(covariates[,c("intercept","X1")]) %*% t(t(beta1))
  mu2 <- as.matrix(covariates[,c("intercept","X2")]) %*% t(t(beta2))
  
  # UforY <- runif(nsample,0,1)
  
  U1 <- unlist(lapply(exp(mu1),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  U2 <- unlist(lapply(exp(mu2),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  Y <- (U1 > 0) *  U2
  
  ### Measurement error model
  
  
  muZplus <- as.matrix(covariates[,c("intercept")]) %*% t(t(alphaplus))
  muZminus <- as.matrix(covariates[,c("intercept")]) %*% t(t(alphaminus))
  # muZplus <- as.matrix(covariates[,c("Xplus")]) * alphaplus
  # muZminus <- as.matrix(covariates[,c("Xminus")]) * alphaminus
  
  Zplus <- unlist(lapply(exp(muZplus),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  
  muZminuspnorm <- expit(muZminus)
  Zminus <- unlist(lapply(1:nsample,FUN=function(x){
    rbinom(1, size = Y[x], prob = muZminuspnorm[x])
  }))
  
  Ystar <- Y + Zplus - Zminus
  
  return(list(main=data.frame(Ystar=Ystar[1:nmain],
                              covariates[1:nmain,],
                              Zplus=Zplus[1:nmain], 
                              Zminus=Zminus[1:nmain]
                              ),
              validation=data.frame(
                    Ystar=Ystar[(nmain+1):(nmain+nvalid)],
                    Y=Y[(nmain+1):(nmain+nvalid)],
                    covariates[(nmain+1):(nmain+nvalid),c("intercept")],
                    Zplus=Zplus[(nmain+1):(nmain+nvalid)],
                    Zminus=Zminus[(nmain+1):(nmain+nvalid)],
                    muZminuspnorm= muZminuspnorm[(nmain+1):(nmain+nvalid)])))
}

### This function is for generating the simulated data with Internal Validation Sample

dataSim5 <- function(seedid, nmain, nvalid, beta1, beta2, alphaplus, alphaminus){
  
  nsample <- nmain + nvalid
  set.seed(seed_i[seedid])
  
  ### Response model 
  
  # X1 <- runif(nsample,0,1)
  X1 <- rbinom(nsample,1,0.5)
  X2 <- runif(nsample,0,1)
  Xplus <-  runif(nsample,0,1)
  Xminus <- runif(nsample,0,2)
  
  covariates <- data.frame(intercept=rep(1,nsample), X1=X1, X2=X2)
  
  beta1 <- beta1
  beta2 <- beta2
  
  mu1 <- as.matrix(covariates[,c("intercept","X1")]) %*% t(t(beta1))
  mu2 <- as.matrix(covariates[,c("intercept","X2")]) %*% t(t(beta2))
  
  # UforY <- runif(nsample,0,1)
  
  U1 <- unlist(lapply(exp(mu1),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  U2 <- unlist(lapply(exp(mu2),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  Y <- (U1 > 0) *  U2
  
  ### Measurement error model
  
  
  muZplus <- as.matrix(covariates[,c("intercept")]) * alphaplus
  muZminus <- as.matrix(covariates[,c("intercept")]) * alphaminus
  
  Zplus <- unlist(lapply(exp(muZplus),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  
  muZminuspnorm <- expit(muZminus)
  Zminus <- unlist(lapply(1:nsample,FUN=function(x){
    rbinom(1, size = Y[x], prob = muZminuspnorm[x])
  }))
  
  Ystar <- Y + Zplus - Zminus
  
  return(list(main=data.frame(Ystar=Ystar[1:nmain],
                              covariates[1:nmain,],
                              Zplus=Zplus[1:nmain], 
                              Zminus=Zminus[1:nmain]
  ),
  validation=data.frame(
    Ystar=Ystar[(nmain+1):(nmain+nvalid)],
    Y=Y[(nmain+1):(nmain+nvalid)],
    covariates[(nmain+1):(nmain+nvalid),],
    Zplus=Zplus[(nmain+1):(nmain+nvalid)],
    Zminus=Zminus[(nmain+1):(nmain+nvalid)],
    muZminuspnorm= muZminuspnorm[(nmain+1):(nmain+nvalid)])))
}

### This function is for generating the simulated data of the THIRD simulation studies with continuous covariates

dataSim6 <- function(seedid, nsample, beta1, beta2, alphaplus, alphaminus){
  
  set.seed(seed_i[seedid])
  
  ### Response model 
  
  # X1 <- runif(nsample,0,1)
  X1 <- rbinom(nsample,1,0.5)
  X2 <- runif(nsample,0,1)
  Xplus <-  runif(nsample,0,1)
  # Xminus <- rbinom(nsample,1,0.5)
  Xminus <- runif(nsample,0,2)
  
  covariates <- data.frame(intercept=rep(1,nsample),X1=X1, X2=X2, Xplus=Xplus, Xminus=Xminus)
  
  beta1 <- beta1
  beta2 <- beta2
  
  mu1 <- as.matrix(covariates[,c("intercept","X1")]) %*% t(t(beta1))
  mu2 <- as.matrix(covariates[,c("intercept","X2")]) %*% t(t(beta2))
  
  # UforY <- runif(nsample,0,1)
  
  U1 <- unlist(lapply(exp(mu1),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  U2 <- unlist(lapply(exp(mu2),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  Y <- (U1 > 0) *  U2
  
  ### Measurement error model
  
  
  muZplus <- as.matrix(covariates[,"intercept"]) * alphaplus
  muZminus <- as.matrix(covariates[,"intercept"]) * alphaminus
  # muZplus <- as.matrix(covariates[,c("Xplus")]) * alphaplus
  # muZminus <- as.matrix(covariates[,c("Xminus")]) * alphaminus
  
  Zplus <- unlist(lapply(exp(muZplus),FUN=function(x){
    rpois(1,lambda = x)
  }))
  
  
  muZminuspnorm <- expit(muZminus)
  Zminus <- unlist(lapply(1:nsample,FUN=function(x){
    rbinom(1, size = Y[x], prob = muZminuspnorm[x])
  }))
  
  Ystar <- Y + Zplus - Zminus
  
  return(data.frame(Ystar=Ystar,Y=Y,covariates,Zplus=Zplus,Zminus=Zminus,U1=U1,U2=U2))
}


##### 3.2 Functions Generating Gibbs Step #####


fY <- function(Yi, phi, mu2){
  if ( Yi == 0 ){
    return(1-phi+phi*(dpois(0,exp(mu2))))
  } else {
    return(phi * (dpois(Yi,exp(mu2))) )
  }
}

fZplusYZminus <- function(Zplus, Yi, Zminus, phi, mu2, muZplus, muZminus){
  fY(Yi,phi,mu2) * dpois(Zplus, exp(muZplus)) * dbinom(Zminus, Yi, pnorm(muZminus))
}


GenerateJoint <- function(Yistar, phii, mu2i, muZplusi, muZminusi){
  u <- runif(1,0,1)
  corPlus <- NULL
  corMinus <- NULL
  P <- NULL
  for (Zplusi in 0:Yistar){
    for (Zminusi in 0:15) {
      corPlus <- c(corPlus,Zplusi)
      corMinus <- c(corMinus,Zminusi)
      P <- c(P, fZplusYZminus(Zplusi, Yistar-Zplusi+Zminusi, Zminusi, phii, mu2i, muZplusi, muZminusi))
    }
  }
  
  index <- min(which(u<=cumsum(P/sum(P))))
  return(c(corPlus[index], corMinus[index], Yistar - corPlus[index] + corMinus[index]))
}



GenerateAlphaplus <- function(alphaplus, covplus, Zplus){
  alphaplusgen<- rep(0,length(alphaplus))
  for (i in 1:length(alphaplus)) {
    covsumtilde <-  as.matrix(covplus[,-i]) %*% t(t(alphaplus[-i])) 
    shape <- sum(Zplus[covplus[,i]==1],na.rm=T) + 1
    rate <-  sum(exp(covsumtilde[covplus[,i]==1]),na.rm=T)
    alphaplusgen[i] <- log(rgamma(1,shape,rate=rate))
  }
  return(alphaplusgen)
}


GeneratePoiPar <- function(Par, covariates, Outcome, priorgamma){
  Pargen <- rep(0,length(Par))
  for (i in 1:length(Par)) {
    covsumtilde <-  as.matrix(covariates[,-i]) %*% t(t(Par[-i])) 
    shape <- sum(Outcome[covariates[,i]==1],na.rm=T)  + priorgamma[1] # + 1
    rate <-  sum(exp(covsumtilde[covariates[,i]==1]),na.rm=T) + priorgamma[2]
    Pargen[i] <- log(rgamma(1,shape,rate=rate))
  }
  return(Pargen)
}

### 3.3 Generate the Parameter for the Binomial models  ####

gTruncRobert <- function(trunca){
  alphastar <- (trunca + sqrt(trunca^2 + 4))/2
  z <- rexp(1,alphastar) + trunca
  rho <- exp(-(z-alphastar)^2/2)
  u <- runif(1,0,1) 
  while(u > rho ){
    z <- rexp(1,alphastar) + trunca
    rho <- exp(-(z-alphastar)^2/2)
    u <- runif(1,0,1) 
  }
  return(z)
}

gTruncNorm <- function(a,mu,sign){
  if (sign>0) {
    if (mu >= a) {
      u <- rnorm(1,mu,1)
      while (u < a ) {
        u <- rnorm(1,mu,1)
      }
      return(u)
    } else{
      Pos <- gTruncRobert(a-mu)
      return(Pos+mu)
    }
  } else{
    Neg <- gTruncNorm(-a,-mu,-sign)
    return(-Neg)
  }
}

# 3.4 Generate Hidden Variable  ####

GenerateV <- function(Yi,Zminusi, muZminusi, covminusi){
  Vcovar <- NULL
  V <- NULL
  if (Yi>0){
    if (Zminusi==0) {
      for (j in 1:Yi){
        Vi <- gTruncNorm(0,muZminusi,-1)
        V <- c(V,Vi)
        Vcovar <- rbind(Vcovar,covminusi)
      }
    } else {
      if (Zminusi == Yi){
        for (j in 1:Yi){
          Vi <- gTruncNorm(0,muZminusi,1)
          V <- c(V,Vi)
          Vcovar <- rbind(Vcovar,covminusi)
        }
      } else {
        for (j in 1:Zminusi){
          Vi <- gTruncNorm(0,muZminusi,1)
          V <- c(V,Vi)
          Vcovar <- rbind(Vcovar,covminusi)
        }
        for (j in (Zminusi+1):Yi){
          Vi <- gTruncNorm(0,muZminusi,-1)
          V <- c(V,Vi)
          Vcovar <- rbind(Vcovar,covminusi)
        }
      }
    }
  }
  return(list(V=V, Vcovar=Vcovar))
}

GenerateAlphaminus <- function(alphaminus, Y, Covarminus, Zminus, nsample, priormu, priorSigma){
  muZminusgen <- as.matrix(Covarminus) %*% t(t(alphaminus)) 
  
  newV <- lapply(1:nsample,FUN = function(i){
    GenerateV(Y[i], Zminus[i], muZminusgen[i], Covarminus[i,])
  })
  
  V <- unlist(lapply(newV, `[`, "V"))
  
  Vcovartemp <- do.call(mapply, c(FUN=rbind, lapply(newV, `[`, "Vcovar")))
  if (dim(Covarminus)[2]==1) {Vcovar <- Vcovartemp} else {Vcovar <- do.call(cbind, Vcovartemp)}
  Vcovar <- setNames(Vcovar,names(Covarminus))
  
  alphaSigma <- solve(solve(priorSigma)+ t(Vcovar) %*% Vcovar)
  alphamu <- alphaSigma %*% ( solve(priorSigma) %*% priormu + t(Vcovar) %*% V)
  alphaminus_gen <- mvrnorm(1,mu = alphamu, Sigma=alphaSigma )
  return(alphaminus_gen)
}


## GenerateAlphaminus c++ version
GenerateAlphaminusCpp <- function(alphaminus, Y, Covarminus, Zminus, priormu, priorSigma){
  muZminusgen <- as.matrix(Covarminus) %*% t(t(alphaminus)) 
  
  Vmat <- ZI_GenerateV(Y, Zminus, muZminusgen,  as.matrix(Covarminus))
  
  V <- Vmat[,1]
  Vcovar <- Vmat[,-1]
  
  # alphaSigma <- solve(t(Vcovar) %*% Vcovar)
  # alphamu <- alphaSigma %*% ( t(Vcovar) %*% V)
  
  alphaSigma <- solve(solve(priorSigma)+ t(Vcovar) %*% Vcovar)
  alphamu <- alphaSigma %*% ( solve(priorSigma) %*% priormu + t(Vcovar) %*% V)
  alphaminus_gen <- mvrnorm(1,mu = alphamu, Sigma=alphaSigma )
  return(alphaminus_gen)
}


gTruncPois <- function(lambda){
  xstar <- rpois(1,lambda=lambda)
  
  while (xstar==0) {
    xstar <- rpois(1,lambda=lambda)
  }
  
  return(xstar)
}



generateZs <- function(yistar, yi, muZminusi){
  Zmcandidate <- rbinom(1,yi,expit(muZminusi))
  Zpcandidate <- yistar - yi + Zmcandidate
  if (Zpcandidate>=0) {return(list(Zm=Zmcandidate, Zp=Zpcandidate))}
   else {generateZs(yistar, yi, muZminusi)}
}


GenerateU1 <- function(Y, U2, mu1, nsample){
  U1_gen <- unlist(lapply(1:nsample, FUN=function(i){
    if (Y[i]>0) {return(gTruncPois(mu1[i]))} else{
      if (U2[i]>0) {return(0)} else {
        return(rpois(1,lambda=mu1[i]))
      }
    }
  }))
  return(U1_gen)
}

GenerateU2 <- function(Y, U1, mu2, nsample){
  U2_gen <- unlist(lapply(1:nsample, FUN=function(i){
    if (U1[i]==0) {return(rpois(1,lambda=mu2[i]))} else return(Y[i])
  }))
  return(U2_gen)
}

GenerateU1i <- function(Yi, U2i, mu1i){
  mu1iexp <- exp(mu1i)
  if (Yi>0) {return(gTruncPois(mu1iexp))} else{
    if (U2i>0) {return(0)} else {
      return(rpois(1,lambda=mu1iexp))
    }
  }
}

GenerateU2i <- function(Yi, U1i, mu2i){
  mu2iexp <- exp(mu2i)
  if (U1i==0) {return(rpois(1,lambda=mu2iexp))} else return(Yi)
}

## Metropolis Method
jumpfunc <- function(theta,sigma){
  return(rnorm(1,theta,sigma))
}

GenerateBetaMetro <- function(beta, U, Covar, propsigma, priorgamma){
  mu <-  as.matrix(Covar) %*% t(t(beta))
  
  for (i in  1:length(beta)){
    parastar <- jumpfunc(beta[i],propsigma[i])
    
   
    muprop <- mu + Covar[,i]*(parastar-beta[i])
     
    
    logratioposterior <- (t(U) %*% as.matrix(Covar[,i])+priorgamma[1]-1) * (parastar - beta[i]) +
      sum(exp(mu)-exp(muprop),na.rm=T) + priorgamma[2] * (exp(beta[i])-exp(parastar))
    
    cat("P=",logratioposterior , " ")
    
    if (logratioposterior>=0) {
      beta[i] <- parastar
      mu <- muprop
    } else{
      u <- runif(1,0,1)
      if (log(u)<=logratioposterior) {
        beta[i] <- parastar
        mu <- muprop
      }
    }
  }
  return(beta)
}


#### 4. Diagnostic Tools ####

RejectionRate <- function(sequence){
  sequence1 <- sequence[-length(sequence)]
  sequence2 <- sequence[-1]
  return( sum(sequence1 == sequence2,na.rm= T) / (length(sequence)-1) )
}

autosigma <- function(sigma, results, benchmark){
  sigmapro <- sigma
  index <- 0
  results <- as.matrix(results)
  cat("Rejection Rates: ")
  for (i in 1:dim(results)[2]){
    AR <- RejectionRate(results[,i])
    cat(AR, " ")
    sigmapro[i] <-  sigma[i] - sigma[i] * (AR - benchmark) * 0.5
  }
  
  return(sigmapro)
}

purifyseq <- function(seq, burnin, itskip){
  return(seq[burnin + 1:floor((length(seq)-burnin)/itskip) *itskip])
}

getHDP <- function(seq, alpha){
  nleng<- length(seq)
  seq.sort <- sort(seq)
  IntCandidate <- NULL
  for (i in 1:(nleng-floor((1-alpha)*nleng))){
    IntCandidate <- c(IntCandidate,seq.sort[i+floor((1-alpha)*nleng)]-seq.sort[i])
  }
  index <- which(IntCandidate==min(IntCandidate,na.rm=T))
  if (length(index)>1) {index <- index[1]}
  return(c(seq.sort[index],seq.sort[index+floor((1-alpha)*nleng)]))
}


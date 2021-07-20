### This file shows the example of running one iteration in the proposed method in Simulation Study 3 (Section 6.3)

#### 1. Global Parameters  ####

seed <- 2020
set.seed(seed)

library(MASS)
library(GeneErrorMis)
library(doParallel)

#### 2. Auxiliary Functions  ####
source("Functions.R")

#### 3. Main Implementation Functions ####

Mainfunction_concpp <- function(Ystar, Covarmain1, Covarmain2, Covarplus, Covarminus,
                             beta1, beta2, alphaplus, alphaminus,
                             Uibound,
                             priorgamma,priormu,priorSigma, propsigmai,
                             seed, nmcmc){
  set.seed(seed)
  library(MASS)
  library(GeneErrorMis)
  source("Functions.R")
  nsample <- length(Ystar)
  
  Zplus <- rep(0,nsample)
  Zminus<- rep(0,nsample)
  
  mu1 <-  as.matrix(Covarmain1) %*% t(t(beta1))
  mu2 <- as.matrix(Covarmain2) %*% t(t(beta2))
  muZplus <- as.matrix(Covarplus) * alphaplus 
  muZminus <- as.matrix(Covarminus) * alphaminus
  U1 <- rep(0,nsample)
  U2 <- rep(0,nsample)
  
  Y <- Ystar - Zplus + Zminus
  
  beta1record <- beta1
  beta2record <- beta2
  alphaplusrecord <- alphaplus
  alphaminusrecord <- alphaminus
  for (repk in 1:nmcmc){
    
    if (repk %% 500 ==0){
      cat(repk," ")
    }
    
    
    ### Step 1: Measurement Error: Positive Process
    
    for (t in 1:nsample){
      NewGen <-  ZI_GenerateBigJoint(Ystar[t], Uibound[1], Uibound[2],  mu1[t], mu2[t], muZplus[t],  muZminus[t])
      Zplus[t] <- NewGen[1]
      Zminus[t] <- NewGen[2]
      Y[t] <- NewGen[3]
      U1[t] <- NewGen[4]
      U2[t] <- NewGen[5]
    }
    
    alphaplussafe <- alphaplus
    
    alphaplus <- ZI_GenerateBetaMetro(alphaplus, Zplus, as.matrix(Covarplus), propsigmai,  priorgamma)
    
    if (any(alphaplus< -5)) {alphaplus <- alphaplussafe}
    
    muZplus <- as.matrix(Covarplus) * alphaplus
    
    
    # ### Step 2: Measurement Error: Negative Process
    
    alphaminussafe <- alphaminus
    
    alphaminus <- ZI_GenerateAlphaMNMetro(alphaminus, Y, Zminus, as.matrix(Covarminus),  propsigma= propsigmai,  priormu,  priorSigmas=priorSigma)
    
    if (any(alphaminus< -5)) {alphaminus <- alphaminussafe}
    muZminus <- as.matrix(Covarminus) * alphaminus
    
    # ### Step 3: Main Model
    
    ### Step 3.1 Update beta1 and beta 2
    
    beta1 <- ZI_GenerateBetaMetro(Par=beta1, U=U1, Covar=as.matrix(Covarmain1), propsigma=propsigmai[1:2], priorgamma)
    mu1 <-  as.matrix(Covarmain1) %*% t(t(beta1))
    
    beta2 <- ZI_GenerateBetaMetro(Par=beta2, U=U2, Covar=as.matrix(Covarmain2), propsigma=propsigmai[3:4], priorgamma)
    
    mu2 <- as.matrix(Covarmain2) %*% t(t(beta2))
    
    
    beta1record <- rbind(beta1record,beta1)
    beta2record <- rbind(beta2record,beta2)
    alphaplusrecord <- c(alphaplusrecord,alphaplus)
    alphaminusrecord <- c(alphaminusrecord,alphaminus)
  }
  
  BayesResults <- data.frame(beta1record,beta2record,alphaplusrecord,alphaminusrecord)

  
  ### Summarize results from sample MCMC data
  
  results <- apply(BayesResults, MARGIN =2, FUN = function(x){
    xpure <- purifyseq(x, 500, 5)
    mean <- mean(xpure,na.rm=T)
    sd <- sd(xpure,na.rm=T)
    median <- median(xpure,na.rm=T)
    CrI <- quantile(xpure, c(0.025,0.975),na.rm=T)
    HDP <- getHDP(xpure,0.05)
    modeH <- getHDP(xpure,0.95)
    return(c(mean,median,sd,CrI,HDP,modeH))
  })

  return(results)
}


#### 4. Simulation Start ####

## var2: Identification of the setting
## var2 = 1: improper prior
## var2 = 0: proper prior

Simbin_gam <- lapply(1:2, FUN=function(var2){
  start <- proc.time()
  
  ## Generate the data
  datasim <- dataSim6(seedid= seed, nsample = 1000, beta1 = c(-0.7,0.7), beta2 = c(1,-0.5), alphaplus = -1, alphaminus = -1)
  
  
  ## Implement the proposed methods
  result <-  Mainfunction_concpp(Ystar = datasim$Ystar,
                                 Covarmain1 = datasim[,c("intercept","X1")],
                                 Covarmain2 = datasim[,c("intercept","X2")],
                                 Covarplus = datasim[,c("intercept")],
                                 Covarminus = datasim[,c("intercept")],
                                 beta1 = c(-0.5,0.5), beta2 = c(0,0),    
                                 alphaplus = c(0), alphaminus=c(0),
                                 Uibound = c(7,11),
                                 priorgamma = ifelse(var2==1,c(0.001^2,1000^2),c(1,1)),
                                 priormu =  -1,
                                 priorSigma = ifelse(var2==1,1000^2,1), propsigmai = c(rep(0.05,4),0.05,0.15),
                                 seed = seed*100, nmcmc = 2000)  
  end <- proc.time() - start
  cat(paste0("task",var2,":",end[3],"\n"))
  return(result)
})

save(Simbin_gam,file=paste0("ExampleResults.RData"))


#### 5 Data Structure Explanation ####

# The output object "Simbin_gam" is a list contains two elements:
#     elements 1 corresponds to the results with scenario of improper prior and elements 2 corresponds to the results of the scenario with proper prior
# Each elements composed a 9 X 6 tables:
# The rows corresponds to the measures:
# "mean", "standerd error", "median", "credible interval(Lower Bound)", "credible interval(Upper Bound)",
# "high density credible interval(Lower Bound)", "high density credible interval(Upper Bound)", and "mode"

# The column corresponds to the parameters:
# "beta_phi0", "beta_phix", "beta_mu0", "\beta_mux", "alpha_+0", "alpha_-0"

Simbin_gam



### This simulation studies study the case with different prior choices
### under the case where the covariates of main model are continuous.

### Simulation 3 is step further from Simulation 2 where we try more degree of measurement error


#### 1. Global Parameters  ####

ncores <- 40
ProjectName <- "zeroInf"
# nsample <- 2000

set.seed(2020)
seed_i <- sample(1000000,1000)


range <- c(1,400)

library(MASS)
library(GeneErrorMis)
library(doParallel)

#### 2. Auxiliary Functions  ####
source("code/Simulation/Functions.R")

#### 3. Main Implementation Functions ####

Mainfunction_concpp <- function(Ystar, Covarmain1, Covarmain2, Covarplus, Covarminus,
                             beta1, beta2, alphaplus, alphaminus,
                             Uibound,
                             priorgamma,priormu,priorSigma, propsigmai,
                             seed, nmcmc){
  set.seed(seed)
  library(MASS)
  library(GeneErrorMis)
  source("code/Simulation/Functions.R")
  nsample <- length(Ystar)
  
  Zplus <- rep(0,nsample)
  Zminus<- rep(0,nsample)
  
  mu1 <-  as.matrix(Covarmain1) %*% t(t(beta1))
  mu2 <- as.matrix(Covarmain2) %*% t(t(beta2))
  muZplus <- as.matrix(Covarplus) %*% t(t(alphaplus)) 
  muZminus <- as.matrix(Covarminus) %*% t(t(alphaminus))
  U1 <- rep(0,nsample)
  U2 <- rep(0,nsample)
  
  Y <- Ystar - Zplus + Zminus
  
  beta1record <- beta1
  beta2record <- beta2
  alphaplusrecord <- alphaplus
  alphaminusrecord <- alphaminus
  for (repk in 1:nmcmc){
    
    ## Adjust the proposal standard error 
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
    
    alphaplus <- ZI_GenerateBetaMetro(alphaplus, Zplus, as.matrix(Covarplus), propsigmai[5:6],  priorgamma)
    
    if (any(alphaplus< -5)) {alphaplus[alphaplus< -5] <- alphaplussafe[alphaplus< -5]}
    
    muZplus <- as.matrix(Covarplus) %*% t(t(alphaplus))
    
    
    # ### Step 2: Measurement Error: Negative Process
    
    alphaminussafe <- alphaminus
    
    alphaminus <- ZI_GenerateAlphaMNMetro(alphaminus, Y, Zminus, as.matrix(Covarminus),  propsigma= propsigmai[7:8],  priormu,  priorSigmas=priorSigma)
    
    if (any(alphaminus< -5)) {alphaminus[alphaminus< -5] <- alphaminussafe[alphaminus< -5]}
    
    muZminus <- as.matrix(Covarminus) %*% t(t(alphaminus))
    
    # ### Step 3: Main Model
    
    ### Step 3.1 Update beta1 and beta 2
    
    beta1 <- ZI_GenerateBetaMetro(Par=beta1, U=U1, Covar=as.matrix(Covarmain1), propsigma=propsigmai[1:2], priorgamma)
    mu1 <-  as.matrix(Covarmain1) %*% t(t(beta1))
    
    beta2 <- ZI_GenerateBetaMetro(Par=beta2, U=U2, Covar=as.matrix(Covarmain2), propsigma=propsigmai[3:4], priorgamma)
    
    mu2 <- as.matrix(Covarmain2) %*% t(t(beta2))
    
    
    beta1record <- rbind(beta1record,beta1)
    beta2record <- rbind(beta2record,beta2)
    alphaplusrecord <- rbind(alphaplusrecord,alphaplus)
    alphaminusrecord <- rbind(alphaminusrecord,alphaminus)
  }
  
  BayesResults <- data.frame(beta1record,beta2record,alphaplusrecord,alphaminusrecord)

 
  results <- apply(BayesResults, MARGIN =2, FUN = function(x){
    xpure <- purifyseq(x, 500, 15)
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

remove_outliers_label <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  ifelse(is.na(y),NA,1)
}




#### 4. Main Implementation Functions ####

datasim <- datagencontinue(seedid= 2, nsample = 5000, beta1 = c(-0.7,0.7), beta2 = c(1,-0.5), alphaplus = c(0.5), alphaminus=-2.3)

Results0 <- Mainfunction_concpp(Ystar = datasim$Ystar,
                            Covarmain1 = datasim[,c("intercept","X1")],
                            Covarmain2 = datasim[,c("intercept","X2")],
                            Covarplus = datasim[,c("Xplus")],
                            Covarminus = datasim[,c("Xminus")],
                             beta1 = c(-0.7,0.7), beta2 = c(1,-0.5),
                            alphaplus = c(0.5), alphaminus=-2.3,
                            Uibound = c(7,11),
                            priorgamma = c(0.001,0.001),
                            priormu = 0,
                            priorSigma = 1000, propsigmai = c(rep(0.05,6)),
                             seed = 1, nmcmc = 20000) #120000

Results <- Results0
Results <-  result
RejectionRate(Results[,1])
RejectionRate(Results[,2])
RejectionRate(Results[,3])
RejectionRate(Results[,4])
RejectionRate(Results[,5])
RejectionRate(Results[,6])
RejectionRate(Results[,7])
RejectionRate(Results[,8])

Results <- apply(Results0, MARGIN =2, FUN = function(x){
  xpure <- purifyseq(x, 500, 50)
  return(xpure)
})

nsample <- dim(Results)[1]
colMeans(Results[1:nsample,1:2])
xcord <- 1:nsample
plot(Results[,1]~xcord,type="l")
quantile(Results[500:dim(Results)[1],1], c(0.025,0.975))
plot(Results[,2]~xcord,type="l")
quantile(Results[500:dim(Results)[1],2], c(0.025,0.975))

# beta2   # c(1.5,-0.5)
colMeans(Results[1:nsample,3:4])
xcord <- 1:nsample
plot(Results[,3]~xcord,type="l")
quantile(Results[500:dim(Results)[1],3], c(0.025,0.975))
plot(Results[,4]~xcord,type="l")
quantile(Results[500:dim(Results)[1],4], c(0.025,0.975))

# alphaplus  # c(-1,0.5)
colMeans(Results[1:nsample,5:7])
xcord <- 1:nsample
plot(Results[xcord,5]~xcord,type="l")
quantile(Results[1:dim(Results)[1],5], c(0.025,0.975))
plot(Results[,6]~xcord,type="l")
quantile(Results[1:dim(Results)[1],6], c(0.025,0.975))
plot(Results[,7]~xcord,type="l")
quantile(Results[1:dim(Results)[1],6], c(0.025,0.975))
# alphaminus  # c(-1)
colMeans(Results[1:nsample,5:6])
xcord <- 1:nsample
plot(Results[xcord,5]~xcord,type="l")
quantile(Results[1:dim(Results)[1],5], c(0.025,0.975))




#### 5. Simulation Start ####

 
alphaminus.candidate <- matrix(c(-1,-1.2,-0.8,-1), ncol = 2)
alphaplus.candidate <- matrix(c(-1,0.6,2,-1.2), ncol = 2)


### 5.1 Study 1 - low addup process  ####

ProjectName <- paste0("ZI5000_Sim3_lowaddup_",range[1],"-",range[2])

cl = makeCluster(ncores, type = "PSOCK",outfile=paste0(ProjectName,".txt"))
clusterExport(cl=cl, varlist=c("seed_i", "nsample"))
registerDoParallel(cl)


Simbin_gam <- foreach(var2=c(1:2)) %:%  foreach(var1=range[1]:range[2]) %dopar% {
  start <- proc.time()
  datasim <- dataSim3(seedid= var1, nsample = 5000, beta1 = c(-0.2,0.7), beta2 = c(1,0.5), alphaplus = alphaplus.candidate[,1], alphaminus=alphaminus.candidate[,var2])
  result <-  Mainfunction_concpp(Ystar = datasim$Ystar,
                                 Covarmain1 = datasim[,c("intercept","X1")],
                                 Covarmain2 = datasim[,c("intercept","X2")],
                                 Covarplus = datasim[,c("intercept","Xplus")],
                                 Covarminus = datasim[,c("intercept","Xminus")],
                                 beta1 = c(0.7,-0.7), beta2 = c(1,-1.5),    # was c(1.2,-1.5) for 1-240
                                 alphaplus = c(0,0), alphaminus=c(0,0),
                                 Uibound = c(7,11),
                                 priorgamma = c(0.001,0.001),
                                 priormu = alphaminus.candidate[,var2],
                                 priorSigma = c(1,1), propsigmai = c(rep(0.05,4),0.05,0.12,0.15,0.2),
                                 seed = var1, nmcmc = 10000)   # , error = function(e) NULL)
  end <- proc.time() - start
  cat(paste0(var1, ": ",end[3],"\n"))
  return(result)
}

save(Simbin_gam,file=paste0("output/Simulation/",ProjectName,".RData"))


### 5.2 Study 2 - high addup process  ####


ProjectName <- paste0("ZI5000_Sim3_highaddup_",range[1],"-",range[2])

cl = makeCluster(ncores, type = "PSOCK",outfile=paste0(ProjectName,".txt"))
clusterExport(cl=cl, varlist=c("seed_i", "nsample"))
registerDoParallel(cl)

Simbin_gam <- foreach(var2=c(1:2)) %:%  foreach(var1=range[1]:range[2]) %dopar% {
  start <- proc.time()
  datasim <- dataSim3(seedid= var1, nsample = 5000, beta1 = c(-0.2,0.7), beta2 = c(1,0.5), alphaplus = alphaplus.candidate[,2], alphaminus=alphaminus.candidate[,var2])
  result <-  Mainfunction_concpp(Ystar = datasim$Ystar,
                                 Covarmain1 = datasim[,c("intercept","X1")],
                                 Covarmain2 = datasim[,c("intercept","X2")],
                                 Covarplus = datasim[,c("intercept","Xplus")],
                                 Covarminus = datasim[,c("intercept","Xminus")],
                                 beta1 = c(0.7,-0.3), beta2 = c(2,-1),
                                 alphaplus = c(0,0), alphaminus=c(1,0),  # was 241-400 alphaminus=c(1,2)
                                 Uibound = c(7,11),
                                 priorgamma = c(0.001,0.001),
                                 priormu = alphaminus.candidate[,var2],
                                 priorSigma = c(1,1), propsigmai = c(rep(0.05,4),0.05,0.12,0.15,0.2),
                                 seed = var1, nmcmc = 10000)   # , error = function(e) NULL)
  end <- proc.time() - start
  cat(paste0(var1, ": ",end[3],"\n"))
  return(result)
}

save(Simbin_gam,file=paste0("output/Simulation/",ProjectName,".RData"))



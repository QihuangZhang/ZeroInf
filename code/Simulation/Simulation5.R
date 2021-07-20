### This simulation studies study the case with internal validation data available

### Simulation 4 is step further from Simulation 3 to consider the case when validation is available 
# (We consider the scenario both X_i and W_i is available)


#### 1. Global Parameters  ####

ncores <- 50
ProjectName <- "zeroInf"

set.seed(2020)
seed_i <- sample(1000000,1000)


range <- c(1,1000)

library(MASS)
library(GeneErrorMis)
library(doParallel)

#### 2. Auxiliary Functions  ####
source("code/Simulation/Functions.R")

#### 3. Main Implementation Functions ####

Mainfunction_intcpp <- function(Ystar, Covarmain1, Covarmain2, Covarplus, Covarminus,
                                Ystarval, Yval, Covarvalmain1, Covarvalmain2, Covarvalplus, Covarvalminus,
                                beta1, beta2, alphaplus, alphaminus,
                                Uibound,
                                priorgamma, priormu, priorSigma, propsigmai,
                                seed, nmcmc){
  set.seed(seed)
  library(MASS)
  library(GeneErrorMis)
  source("code/Simulation/Functions.R")
  nmain <- length(Ystar)
  nval <- length(Ystarval)
  
  nsample <- nmain + nval
  
  Zplus <- rep(0,nsample)
  Zminus<- rep(0,nsample)
  
  Covarmain1all <- rbind(Covarmain1, Covarvalmain1)
  Covarmain2all <- rbind(Covarmain2, Covarvalmain2)
  
  mu1 <-  as.matrix(rbind(Covarmain1all)) %*% t(t(beta1))
  mu2 <- as.matrix(rbind(Covarmain2all)) %*% t(t(beta2))
  
  Covarplusall <- c(Covarplus, Covarvalplus)
  Covarminusall <- c(Covarminus, Covarvalminus)
  
  muZplus <- as.matrix(Covarplusall) * alphaplus 
  muZminus <- as.matrix(Covarminusall) * alphaminus
  U1 <- rep(0,nsample)
  U2 <- rep(0,nsample)
  
  Y <- Ystar - Zplus[1:nmain] + Zminus[1:nmain]
  Y <- c(Y,Yval)
  
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
    
    for (t in 1:nmain){
      NewGen <-  ZI_GenerateBigJoint(Ystar[t], Uibound[1], Uibound[2],  mu1[t], mu2[t], muZplus[t],  muZminus[t])
      Zplus[t] <- NewGen[1]
      Zminus[t] <- NewGen[2]
      Y[t] <- NewGen[3]
      U1[t] <- NewGen[4]
      U2[t] <- NewGen[5]
    }

    for (t in 1:nval){
      NewGen2 <- ZI_GenerateZpZmJoint(Ystarval[t], Yval[t], muZplus[t+nmain], muZminus[t+nmain])
      Zminus[t+nmain] <- NewGen2[1]
      Zplus[t+nmain] <- NewGen2[2]
      
      U1[t+nmain] <- GenerateU1i(Yval[t], U2[t+nmain],  mu1[t+nmain])
      U2[t+nmain] <- GenerateU2i(Yval[t], U1[t+nmain],  mu2[t+nmain])
    }
    
    alphaplussafe <- alphaplus
    
    alphaplus <- ZI_GenerateBetaMetro(alphaplus, Zplus, as.matrix(Covarplusall), propsigmai,  priorgamma)
    
    if (any(alphaplus< -5)) {alphaplus[alphaplus< -5] <- alphaplussafe[alphaplus< -5]}
    
    muZplus <- as.matrix(Covarplusall) %*% t(t(alphaplus))
    
    
    # ### Step 2: Measurement Error: Negative Process
    
    alphaminussafe <- alphaminus
    
    alphaminus <- ZI_GenerateAlphaMNMetro(alphaminus, Y, Zminus, as.matrix(Covarminusall),  propsigma= propsigmai,  priormu,  priorSigmas=priorSigma)
    
    if (any(alphaminus< -5)) {alphaminus[alphaminus< -5] <- alphaminussafe[alphaminus< -5]}
    
    muZminus <- as.matrix(Covarminusall) %*% t(t(alphaminus))
    
    # ### Step 3: Main Model

    ### Step 3.1 Update beta1 and beta 2

    beta1 <- ZI_GenerateBetaMetro(Par=beta1, U=U1, Covar=as.matrix(Covarmain1all), propsigma=propsigmai, priorgamma)
    mu1 <-  as.matrix(Covarmain1all) %*% t(t(beta1))

    beta2 <- ZI_GenerateBetaMetro(Par=beta2, U=U2, Covar=as.matrix(Covarmain2all), propsigma=propsigmai, priorgamma)

    mu2 <- as.matrix(Covarmain2all) %*% t(t(beta2))


    beta1record <- rbind(beta1record,beta1)
    beta2record <- rbind(beta2record,beta2)
    alphaplusrecord <- c(alphaplusrecord,alphaplus)
    alphaminusrecord <- c(alphaminusrecord,alphaminus)
  }
  
  BayesResults <- data.frame(beta1record,beta2record,alphaplusrecord,alphaminusrecord)
  
  
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

remove_outliers_label <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  ifelse(is.na(y),NA,1)
}

#### 4. Main Implementation Functions ####

alphaminus.candidate <- matrix(c(-1,-1.2,-0.8,-1), ncol = 2)
alphaplus.candidate <- matrix(c(-1,0.6,2,-1.2), ncol = 2)

datasim <- dataSim5(seedid= 1, nmain = 2500, nvalid = 2500, beta1 = c(-0.2,0.7), beta2 = c(1,0.5), alphaplus = alphaplus.candidate[,1], alphaminus=alphaminus.candidate[,1])

Results0 <- Mainfunction_intcpp(Ystar = datasim$main$Ystar,
                                Covarmain1 = datasim$main[,c("intercept","X1")],
                                Covarmain2 = datasim$main[,c("intercept","X2")],
                                Covarplus = datasim$main[,c("intercept","Xplus")],
                                Covarminus = datasim$main[,c("intercept","Xminus")],
                                Ystarval = datasim$val$Ystar,
                                Yval = datasim$val$Y,
                                Covarvalmain1 = datasim$val[,c("intercept","X1")],
                                Covarvalmain2 = datasim$val[,c("intercept","X2")],
                                Covarvalplus = datasim$val[,c("intercept","Xplus")],
                                Covarvalminus = datasim$val[,c("intercept","Xminus")],
                                beta1 = c(-0.2,0.7), beta2 = c(1,0.5),
                                alphaplus = alphaplus.candidate[,1], alphaminus=alphaminus.candidate[,1],
                                Uibound = c(7,11),
                                priorgamma = c(0.001,0.001),
                                priormu = c(0,0),
                                priorSigma = c(1000,1000), propsigmai = c(c(0.05,0.06,0.03,0.04,0.05,0.1,0.08,0.1)),
                                seed = 1, nmcmc = 5000)

Results <- Results0

RejectionRate(Results[,1])
RejectionRate(Results[,2])
RejectionRate(Results[,3])
RejectionRate(Results[,4])
RejectionRate(Results[,5])
RejectionRate(Results[,6])
RejectionRate(Results[,7])
RejectionRate(Results[,8])


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
xcord <- 1:length(Results[,5])
plot(Results[,5]~xcord,type="l")
quantile(Results[1:dim(Results)[1],5], c(0.025,0.975))
plot(Results[,6]~xcord,type="l")
quantile(Results[1:dim(Results)[1],6], c(0.025,0.975))
plot(Results[,7]~xcord,type="l")
quantile(Results[1:dim(Results)[1],7], c(0.025,0.975))

# alphaminus  # c(-1)
colMeans(Results[1:nsample,5:6])
xcord <- 1:nsample
plot(Results[,8]~xcord,type="l")
quantile(Results[1:dim(Results)[1],8], c(0.025,0.975))




#### 5. Simulation Start ####

# ### Jul 1



### 5.1 Study 1 - low addup process  ####

ProjectName <- paste0("ZI2000_Sim5_",range[1],"-",range[2])

cl = makeCluster(ncores, type = "PSOCK",outfile=paste0(ProjectName,".txt"))
clusterExport(cl=cl, varlist=c("seed_i"))
registerDoParallel(cl)


Simbin_gam <- foreach(var2=c(1:2)) %:%  foreach(var1=range[1]:range[2]) %dopar% {
  start <- proc.time()
  datasim <- dataSim5(seedid= var1, nmain = 1000, nvalid = 1000, beta1 = c(-0.7,0.7), beta2 = c(1,-0.5), alphaplus = -1, alphaminus= -1)
  result <-  tryCatch(Mainfunction_intcpp(Ystar = datasim$main$Ystar,
                                 Covarmain1 = datasim$main[,c("intercept","X1")],
                                 Covarmain2 = datasim$main[,c("intercept","X2")],
                                 Covarplus = datasim$main[,c("intercept")],
                                 Covarminus = datasim$main[,c("intercept")],
                                 Ystarval = datasim$val$Ystar, 
                                 Yval = datasim$val$Y,
                                 Covarvalmain1 = datasim$val[,c("intercept","X1")],
                                 Covarvalmain2 = datasim$val[,c("intercept","X2")],
                                 Covarvalplus = datasim$val[,c("intercept")],
                                 Covarvalminus = datasim$val[,c("intercept")],
                                 beta1 =  c(-0.5,0.5), beta2 = c(0.5,0),
                                 alphaplus = 0, alphaminus= 0,
                                 Uibound = c(7,11),
                                 priorgamma = ifelse(var2==1,c(0.001,0.001),c(1,1)),
                                 priormu = -1,
                                 priorSigma = ifelse(var2==1,1000,1), propsigmai = c(c(0.05,0.06,0.05,0.05,0.05,0.08)),
                                 seed = 1, nmcmc = 6000), error= function(e) return(NA))
  end <- proc.time() - start
  cat(paste0(var1, ": ",end[3],"\n"))
  return(result)
}

save(Simbin_gam,file=paste0("output/Simulation/",ProjectName,".RData"))
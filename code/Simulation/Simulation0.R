### This file is run for the naive model.


#### 1. Global Parameters  ####

ncores <- 40

nsample <- 5000

set.seed(2020)
seed_i <- sample(1000000,1000)

range <- c(1,500)



library(doParallel)

#### 2. Auxiliary Functions  ####
source("code/Simulation/Functions.R")

#### 3. Main Implementation Functions ####
Mainbinnaive <- function(Ystar, Covarmain1, Covarmain2, Covarplus, Covarminus,
                                             beta1, beta2, 
                                             priorgamma, propsigmai,
                                             seed, nmcmc){
  set.seed(seed)
  library(MASS)
  library(GeneErrorMis)
  source("code/Simulation/Functions.R")
  
  
  nsample <- length(Ystar)
  
  
  mu1 <-  as.matrix(Covarmain1) %*% t(t(beta1))
  phi <- cexpexp(as.matrix(Covarmain1) %*% t(t(beta1)))
  mu2 <- as.matrix(Covarmain2) %*% t(t(beta2))
  U2 <- rep(1,length(Ystar)) 
  U1 <- rep(1,length(Ystar)) 
  
  Y <- Ystar
  
  beta1record <- beta1
  beta2record <- beta2
  
  for (repk in 1:nmcmc){
    if (repk %% 500 ==0){
      cat(repk," ")
    }
    ### Step 1: Measurement Error: Positive Process
       
    # ### Step 3: Main Model
    
    ### Step 3.1 Data Augmentation
    U1 <- ZI_GenerateU1(Y,U2,exp(mu1))
    U2 <- ZI_GenerateU2(Y,U1,exp(mu2))
    
    ### Step 3.2 Update beta1 and beta 2
    
    beta1 <- ZI_GenerateBetaMetro(Par=beta1, U=U1, Covar=as.matrix(Covarmain1), propsigma=propsigmai[1:2], priorgamma)
    mu1 <-  as.matrix(Covarmain1) %*% t(t(beta1))
    phi <- cexpexp(mu1)
    
    
    beta2 <- ZI_GenerateBetaMetro(Par=beta2, U=U2, Covar=as.matrix(Covarmain2), propsigma=propsigmai[3:4], priorgamma)
    mu2 <- as.matrix(Covarmain2) %*% t(t(beta2))
    
    
    beta1record <- rbind(beta1record,beta1)
    beta2record <- rbind(beta2record,beta2)
  }
  
  BayesResults <- data.frame(beta1record,beta2record)
  
  results <- apply(BayesResults, MARGIN =2, FUN = function(x){
    xpure <- purifyseq(x, 500, 1)
    mean <- mean(xpure,na.rm=T)
    sd <- sd(xpure,na.rm=T)
    median <- median(xpure,na.rm=T)
    CrI <- quantile(xpure, c(0.025,0.975))
    return(c(mean,median,sd,CrI))
  })
  
  return(results)
}


#### 5. Parallel Computing ####

#### 5.1 Naive Analysis in Simulation 2 Setting 1 (Sample size:5000) ####

ProjectName <- paste0("ZI5000_Naive_Sim2_para1_",range[1],"-",range[2])

cl = makeCluster(ncores, type = "PSOCK",outfile=paste0(ProjectName,".txt"))
clusterExport(cl=cl, varlist=c("seed_i", "nsample"))
registerDoParallel(cl)

# var2 <- 0.001

Simbin_gam0 <- foreach(var2=c(0.001,1)) %:%  foreach(var1=range[1]:range[2]) %dopar% {
  start <- proc.time()
  datasim <- datagencontinue(seedid= var1, nsample = 5000, beta1 = c(-0.7,0.7), beta2 = c(1,-0.5), alphaplus = c(0.5), alphaminus=-2.3)
  result <- tryCatch(Mainbinnaive(seed = var1, Ystar = datasim$Ystar,
                                  Covarmain1 = datasim[,c("intercept","X1")],
                                  Covarmain2 = datasim[,c("intercept","X2")],
                                  Covarplus = datasim[,c("Xplus")],
                                  Covarminus = datasim[,c("Xminus")],
                                  beta1 = c(-0.7,0.7), beta2 = c(1,-0.5),
                                  priorgamma = c(var2,var2),
                                  propsigmai = c(rep(0.05,4)),
                                  nmcmc = 10000), error = function(e) NULL)
  end <- proc.time() - start
  cat(paste0(var1, ": ",end[3],"\n"))
  return(result)
}

save(Simbin_gam0,file=paste0("output/Simulation/",ProjectName,".RData"))


# # #### 5.2 Naive Analysis in Simulation 2 Setting 2 (Sample size:5000) ####
#
ProjectName <- paste0("ZI5000_Naive_Sim2_para2_",range[1],"-",range[2])

cl = makeCluster(ncores, type = "PSOCK",outfile=paste0(ProjectName,".txt"))
clusterExport(cl=cl, varlist=c("seed_i", "nsample"))
registerDoParallel(cl)

# var2 <- 0.001

Simbin_gam0 <- foreach(var2=c(0.001,1)) %:%  foreach(var1=range[1]:range[2]) %dopar% {
  start <- proc.time()
  datasim <- datagencontinue(seedid= var1, nsample = 5000, beta1 = c(-0.2,0.7), beta2 = c(1,0.5), alphaplus = c(0.5), alphaminus=-2.3)
  result <- tryCatch(Mainbinnaive(seed = var1, Ystar = datasim$Ystar,
                                  Covarmain1 = datasim[,c("intercept","X1")],
                                  Covarmain2 = datasim[,c("intercept","X2")],
                                  Covarplus = datasim[,c("Xplus")],
                                  Covarminus = datasim[,c("Xminus")],
                                  beta1 = c(-0.2,0.7), beta2 = c(1,0.5),
                                  priorgamma = c(var2,var2),
                                  propsigmai = c(rep(0.05,4)),
                                  nmcmc = 10000), error = function(e) NULL)
  end <- proc.time() - start
  cat(paste0(var1, ": ",end[3],"\n"))
  return(result)
}


save(Simbin_gam0,file=paste0("output/Simulation/",ProjectName,".RData"))



# #### 5.3 Naive Analysis in Simulation 2 Setting 3 (Sample size:500) ####

#
ProjectName <- paste0("ZI500small_Naive_Sim2_para2_",range[1],"-",range[2])

cl = makeCluster(ncores, type = "PSOCK",outfile=paste0(ProjectName,".txt"))
clusterExport(cl=cl, varlist=c("seed_i", "nsample"))
registerDoParallel(cl)

# var2 <- 0.001

Simbin_gam0 <- foreach(var2=c(1)) %:%  foreach(var1=range[1]:range[2]) %dopar% {
  start <- proc.time()
  datasim <- datagencontinue(seedid= var1, nsample = 500, beta1 = c(-0.2,0.7), beta2 = c(1,0.5), alphaplus = c(0.5), alphaminus=-2.3)
  result <- tryCatch(Mainbinnaive(seed = var1, Ystar = datasim$Ystar,
                                  Covarmain1 = datasim[,c("intercept","X1")],
                                  Covarmain2 = datasim[,c("intercept","X2")],
                                  Covarplus = datasim[,c("Xplus")],
                                  Covarminus = datasim[,c("Xminus")],
                                  beta1 = c(-0.2,0.7), beta2 = c(1,0.5),
                                  priorgamma = c(var2,var2),
                                  propsigmai = c(rep(0.05,4)),
                                  nmcmc = 10000), error = function(e) NULL)
  end <- proc.time() - start
  cat(paste0(var1, ": ",end[3],"\n"))
  return(result)
}


save(Simbin_gam0,file=paste0("output/Simulation/",ProjectName,".RData"))



# #### 5.4 Naive Analysis in Simulation 3 Setting 1 (Sample size:5000) ####

alphaminus.candidate <- matrix(c(-1,-1.2,-0.8,-1), ncol = 2)
alphaplus.candidate <- matrix(c(-1,0.6,2,-1.2), ncol = 2)




ProjectName <- paste0("ZI5000_Naive_Sim3_set1_",range[1],"-",range[2])

cl = makeCluster(ncores, type = "PSOCK",outfile=paste0(ProjectName,".txt"))
clusterExport(cl=cl, varlist=c("seed_i", "nsample"))
registerDoParallel(cl)

# var2 <- 0.001

Simbin_gam0 <- foreach(var2=c(1,2)) %:%  foreach(var1=range[1]:range[2]) %dopar% {
  start <- proc.time()
  datasim <- dataSim3(seedid= var1, nsample = 5000, beta1 = c(-0.2,0.7), beta2 = c(1,0.5), alphaplus = alphaplus.candidate[,1], alphaminus=alphaminus.candidate[,var2])
  result <- tryCatch(Mainbinnaive(seed = var1, Ystar = datasim$Ystar,
                                  Covarmain1 = datasim[,c("intercept","X1")],
                                  Covarmain2 = datasim[,c("intercept","X2")],
                                  Covarplus = datasim[,c("Xplus")],
                                  Covarminus = datasim[,c("Xminus")],
                                  beta1 = c(-0.7,0.7), beta2 = c(1,-0.5),
                                  priorgamma = c(var2,var2),
                                  propsigmai = c(rep(0.05,4)),
                                  nmcmc = 10000), error = function(e) NULL)
  end <- proc.time() - start
  cat(paste0(var1, ": ",end[3],"\n"))
  return(result)
}

save(Simbin_gam0,file=paste0("output/Simulation/",ProjectName,".RData"))



# #### 5.5 Naive Analysis in Simulation 3 Setting 2 (Sample size:5000) ####

ProjectName <- paste0("ZI5000_Naive_Sim3_set2_",range[1],"-",range[2])

cl = makeCluster(ncores, type = "PSOCK",outfile=paste0(ProjectName,".txt"))
clusterExport(cl=cl, varlist=c("seed_i", "nsample"))
registerDoParallel(cl)

# var2 <- 0.001

Simbin_gam0 <- foreach(var2=c(1,2)) %:%  foreach(var1=range[1]:range[2]) %dopar% {
  start <- proc.time()
  datasim <- dataSim3(seedid= var1, nsample = 5000, beta1 = c(-0.2,0.7), beta2 = c(1,0.5), alphaplus = alphaplus.candidate[,2], alphaminus=alphaminus.candidate[,var2])
  result <- tryCatch(Mainbinnaive(seed = var1, Ystar = datasim$Ystar,
                                  Covarmain1 = datasim[,c("intercept","X1")],
                                  Covarmain2 = datasim[,c("intercept","X2")],
                                  Covarplus = datasim[,c("Xplus")],
                                  Covarminus = datasim[,c("Xminus")],
                                  beta1 = c(-0.2,0.7), beta2 = c(1,0.5),
                                  priorgamma = c(var2,var2),
                                  propsigmai = c(rep(0.05,4)),
                                  nmcmc = 10000), error = function(e) NULL)
  end <- proc.time() - start
  cat(paste0(var1, ": ",end[3],"\n"))
  return(result)
}


save(Simbin_gam0,file=paste0("output/Simulation/",ProjectName,".RData"))


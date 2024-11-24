#Study of the benefits of using extra domains
library(ggplot2)
source('utils.R')

doms <- 2:6  #Number of domains considered
N <- 100     #Number of samples generated for each sample size

#Mean and CI for the absolute difference wrt the causal effect
means <- rep(0,length(doms))
confidence <- matrix(0,2,length(doms))

#Parameters
k_U <- 2  #Number of categories of U
k_W <- 2  #Number of categories of W

doX <- 1 #Calculate the interventional distribution do(X=doX)

set.seed(123)

#Start the generation with the first domain
k_E <- 1
PUE <- mat_gen(k_U,k_E)                       #P(U|E)
PWU <- mat_gen(k_W,k_U)                       #P(W|U) 
QU <- vec_gen(k_U,1)                          #Q(U)
QW <- PWU%*%QU                                #Q(W)
PX1U <- runif(k_U)                            #P(X=1|U)
PXU <- array(c(1-PX1U,PX1U),c(1,k_U,2))       #P(X|U)
PxU <- PXU[,,doX+1]                           #P(X=doX|U)
PyUWX <- array(runif(2*k_U*k_W),c(k_W,k_U,2)) #P(Y=1|U,W,X)
PyUWx <- PyUWX[,,doX+1]                       #P(Y=1|U,W,X=doX)
PyUx <- diag(PyUWx%*%PWU)                     #P(Y=1|U,X=doX)

for(l in 1:length(doms)){
  k_E_prev <- k_E
  k_E <- doms[l]
  n <- round(10000/k_E)
  U <- matrix(0,n,k_E)
  W <- matrix(0,n,k_E)
  X <- matrix(0,n,k_E)
  Y <- matrix(0,n,k_E)
  do <- rep(0,N)
  
  PyEx <- rep(0,k_E)
  PWEx <- matrix(0,k_W,k_E)
  
  #Update the P(U|E) matrix by adding new columns (more domains)
  PUE <- cbind(PUE,mat_gen(k_U,k_E-k_E_prev))   #P(U|E)
  
  #Exact values
  PXE <- as.vector(PxU%*%PUE)                   #P(X=doX|E)
  PUEX <- t(t(PxU*PUE)/PXE)                     #P(U|E,X=doX)
  PWEx_t <- PWU%*%PUEX                          #P(W|E,X=doX)
  PyEx_t <- PyUx%*%PUEX                         #P(Y=1|E,X=doX)
  effect <- PyEx_t%*%pseudosolve(PWEx_t)%*%QW   #Q(Y=1|do(X=doX))
  
  for(j in 1:N){
    #Generation of a sample of W in the target domain
    Wnew <- sample(1:k_W,n,replace=TRUE,prob=QW)
    QW_s <- as.vector(table(factor(Wnew, levels=1:k_W)))/n
    
    for(i in 1:k_E){
      U[,i] <- sample(1:k_U,n,replace=TRUE,prob=PUE[,i])
      W[,i] <- sapply(U[,i],substitution)
      X[,i] <- (runif(n)<PxU[U[,i]])*1
      Y[,i] <- (runif(n)<PyUWX[cbind(U[,i],W[,i],X[,i]+1)])*1
      
      #Estimation of P(Y=1|E,X=doX) and P(W|E,X=doX)
      PyEx[i] <- mean(Y[X[,i]==doX,i])
      PWEx[,i] <- as.vector(table(factor(W[X[,i]==doX,i], levels=1:k_W)))/sum(X[,i]==doX)
    }
    #Estimation of the causal effect
    estim_effect <- PyEx%*%pseudosolve(PWEx)%*%QW_s        #Q(Y=1|do(X=1))
    
    #Difference between the estimated and the real values
    do[j] <- estim_effect-effect
  }
  means[l] <- mean(do)
  confidence[1,l] <- quantile(do,0.025)
  confidence[2,l] <- quantile(do,0.975)
}

#Representation of the estimation of the causal effect depending on the sample size
data <- data.frame(doms=doms,means=means,conf1=confidence[1,],conf2=confidence[2,])
data2 <- data[2:nrow(data),]

ggplot(data=data,aes(x=doms)) +
  geom_point(aes(y=means)) +
  geom_errorbar(aes(x=doms,ymin=conf1,ymax=conf2),width=0.2) +
  theme_bw() +
  labs(x='Number of domains', y='Estimation error', color='')


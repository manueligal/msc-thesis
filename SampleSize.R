#Estimation of the causal effect depending on the sample size
library(ggplot2)

#Function to generate k numbers that add up to 1 (2 methods)
vec_gen <- function(k,method){
  if(method==1){
    diff(c(0,sort(runif(k-1)),1))
  }else{
    x <- runif(k)
    x/sum(x)
  }
}

#Function to create a matrix of size kxm with columns that add up to 1
mat_gen <- function(k,m){
  replicate(m,vec_gen(k,1))
}

#Function to obtain a sample of W from a sample of U
substitution <- function(i){
  sample(1:k_W,1,prob=PWU[,i])
}

sample_sizes <- c(200,1000,seq(2e3,1e4,2e3))
N <- 100    #Number of samples generated for each sample size

#Mean and CI for the absolute difference wrt the causal effect
means <- rep(0,length(sample_sizes))
confidence <- matrix(0,2,length(sample_sizes))

#Parameters
k_E <- 2  #Number of domains
k_U <- 2  #Number of categories of U
k_W <- 2  #Number of categories of W

#Generation of the probability matrices
PUE <- mat_gen(k_U,k_E)           #P(U|E)
PWU <- mat_gen(k_W,k_U)           #P(W|U) 
QU <- vec_gen(k_U,1)              #Q(U)
QW <- PWU%*%QU                    #Q(W)

#Exact values
PX1U <- vec_gen(k_U,1)                        #P(X=1|U)
PX1E <- as.vector(PX1U%*%PUE)                 #P(X=1|E)
PUEX1 <- t(t(PX1U*PUE)/PX1E)                  #P(U|E,X=1)
PWEx_t <- PWU%*%PUEX1                         #P(W|E,X=1)
PyUWX <- array(runif(2*k_U*k_W),c(k_W,k_U,2)) #P(Y=1|U,W,X)
PyUWx <- PyUWX[,,2]                           #P(Y=1|U,W,X=1)
PyUx <- diag(PyUWx%*%PWU)                     #P(Y=1|U,X=1)
PyEx_t <- PyUx%*%PUEX1                        #P(Y=1|E,X=1)
effect <- PyEx_t%*%solve(PWEx_t)%*%QW         #Q(Y=1|do(X=1))

for(l in 1:length(sample_sizes)){
  n <- sample_sizes[l]
  U <- matrix(0,n,k_E)
  W <- matrix(0,n,k_E)
  X <- matrix(0,n,k_E)
  Y <- matrix(0,n,k_E)
  do <- rep(0,N)
  
  PyEx <- rep(0,k_E)
  PWEx <- matrix(0,k_W,k_E)
  
  set.seed(123)
  
  for(j in 1:N){
    #Generation of a sample of W in the target domain
    Wnew <- sample(1:k_W,n,replace=TRUE,prob=QW)
    QW_s <- as.vector(table(factor(Wnew, levels=1:k_W)))/n
    
    for(i in 1:k_E){
      U[,i] <- sample(1:k_U,n,replace=TRUE,prob=PUE[,i])
      W[,i] <- sapply(U[,i],substitution)
      X[,i] <- (runif(n)<PX1E[U[,i]])*1
      Y[,i] <- (runif(n)<PyUWX[cbind(U[,i],W[,i],X[,i]+1)])*1
      
      #Estimation of P(Y=1|E,X=1) and P(W|E,X=1)
      PyEx[i] <- mean(Y[X[,i]==1,i])
      PWEx[,i] <- as.vector(table(factor(W[X[,i]==1,i], levels=1:k_W)))/sum(X[,i]==1)
    }
    #Estimation of the causal effect
    estim_effect <- PyEx%*%solve(PWEx)%*%QW_s
    
    #Difference between the estimated and the real values
    do[j] <- estim_effect-effect
  }
  means[l] <- mean(do)
  confidence[1,l] <- quantile(do,0.025)
  confidence[2,l] <- quantile(do,0.975)
}

#Representation of the estimation of the causal effect depending on the sample size
data <- data.frame(size=sample_sizes,means=means,conf1=confidence[1,],conf2=confidence[2,])

ggplot(data=data,aes(x=size)) +
  geom_point(aes(y=means)) +
  geom_errorbar(aes(x=size,ymin=conf1,ymax=conf2),width=150) +
  theme_bw() +
  labs(x='Sample size', y='Estimation error')

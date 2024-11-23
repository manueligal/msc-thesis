#Estimation of the causal effect depending on empirical condition number
library(ggplot2)
source('utils.R')

n <- 5000   #Sample size
M <- 5      #Number of sets of matrices generated
N <- 100    #Number of samples generated for each set of matrices

#Parameters
k_E <- 2  #Number of domains
k_U <- 2  #Number of categories of U
k_W <- 2  #Number of categories of W

doX <- 1 #Calculate the interventional distribution do(X=doX)

data <- data.frame(matrix(NA,M*N,3))
colnames(data) <- c('kappa','do','set')

set.seed(4)

for(m in 1:M){
  #Generation of the probability matrices
  PUE <- mat_gen(k_U,k_E)           #P(U|E)
  PWU <- mat_gen(k_W,k_U)           #P(W|U) 
  QU <- vec_gen(k_U,1)              #Q(U)
  QW <- PWU%*%QU                    #Q(W)
  
  #Exact values
  PX1U <- runif(2)                              #P(X=1|U)
  PXU <- array(c(1-PX1U,PX1U),c(1,2,2))         #P(X|U)
  PxU <- PXU[,,doX+1]                           #P(X=doX|U)
  PXE <- as.vector(PxU%*%PUE)                   #P(X=doX|E)
  PUEX <- t(t(PxU*PUE)/PXE)                     #P(U|E,X=doX)
  PWEx_t <- PWU%*%PUEX                          #P(W|E,X=doX)
  PyUWX <- array(runif(2*k_U*k_W),c(k_W,k_U,2)) #P(Y=1|U,W,X)
  PyUWx <- PyUWX[,,doX+1]                       #P(Y=1|U,W,X=doX)
  PyUx <- diag(PyUWx%*%PWU)                     #P(Y=1|U,X=doX)
  PyEx_t <- PyUx%*%PUEX                         #P(Y=1|E,X=doX)
  effect <- PyEx_t%*%pseudosolve(PWEx_t)%*%QW   #Q(Y=1|do(X=doX))
  
  U <- matrix(0,n,k_E)
  W <- matrix(0,n,k_E)
  X <- matrix(0,n,k_E)
  Y <- matrix(0,n,k_E)
  do <- rep(0,N)
  kappas <- rep(0,N)
  
  PyEx <- rep(0,k_E)
  PWEx <- matrix(0,k_W,k_E)
  
  data$set[(N*(m-1)+1):(N*m)] <- kappa(PWEx_t)
  
  for(j in 1:N){
    #Generation of a sample of W in the target domain
    Wnew <- sample(1:k_W,n,replace=TRUE,prob=QW)
    QW_s <- as.vector(table(factor(Wnew, levels=1:k_W)))/n
    
    for(i in 1:k_E){
      U[,i] <- sample(1:k_U,n,replace=TRUE,prob=PUE[,i])
      W[,i] <- sapply(U[,i],substitution)
      X[,i] <- (runif(n)<PXE[U[,i]])*1
      Y[,i] <- (runif(n)<PyUWX[cbind(U[,i],W[,i],X[,i]+1)])*1
      
      #Estimation of P(Y=1|E,X=doX) and P(W|E,X=doX)
      PyEx[i] <- mean(Y[X[,i]==doX,i])
      PWEx[,i] <- as.vector(table(factor(W[X[,i]==doX,i], levels=1:k_W)))/sum(X[,i]==doX)
    }
    
    #Estimation of the causal effect
    estim_effect <- PyEx%*%pseudosolve(PWEx)%*%QW_s
    
    #Difference between the estimated and the real values
    data$do[N*(m-1)+j] <- estim_effect-effect
    
    #Condition number of the estimation of P(W|E,X=doX)
    data$kappa[N*(m-1)+j] <- kappa(PWEx)
  }
}

ggplot(data=data,aes(x=log(kappa),y=do,color=log(set))) +
  geom_point() +
  labs(x=expression('log('~kappa~'('~widehat(P)~'(W|E,x)))'), y='Estimation error',
       color=expression('log('~kappa~'('~'P(W|E,x)))  ')) +
  theme_bw() +
  theme(legend.direction='vertical',legend.position='right',
        legend.title=element_text(margin=margin(0,0,0,-30,'pt')),
        legend.box.margin=margin(0,0,0,25))

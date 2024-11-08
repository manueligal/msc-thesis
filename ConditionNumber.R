#Estimation of the causal effect depending on the condition number
library(ggplot2)

M <- 20      #Number of different random matrices generated

#Mean and confidence interval for the causal effect
means <- rep(0,M)
confidence <- matrix(0,2,M)
kappas <- rep(0,M)

for(m in 1:M){
  print(m)
  set.seed(100+m)
  p <- runif(2)
  PUE <- cbind(c(p[1],1-p[1]),c(p[2],1-p[2]))/10           #P(U|E)
  PWU <- cbind(c(8,2),c(1,9))/10                           #P(W|U) 
  QW <- c(0.4,0.6)                                         #Q(W)
  
  #Exact values
  PX1U <-  c(6,8)/10                       #P(X=1|U)
  PX1E <- as.vector(PX1U%*%PUE)            #P(X=1|E)
  PUEX1 <- t(t(PX1U*PUE)/PX1E)             #P(U|E,X=1)
  PWEx_t <- PWU%*%PUEX1                    #P(W|E,X=1)
  PyWUx <- rbind(c(0.66,0.68),c(0.7,0.72)) #P(Y=1|W,U,X=1)
  PyUx <- diag(PyWUx%*%PWU)                #P(Y=1|U,X=1)
  PyEx_t <- PyUx%*%PUEX1                   #P(Y=1|E,X=1)
  effect <- PyEx_t%*%solve(PWEx_t)%*%QW    #Q(Y=1|do(X=1))
  
  #Function to obtain a sample of W from a sample of U
  substitution <- function(i){
    sample(1:2,1,prob=PWU[,i])
  }
  
  n <- 1e3
  N <- 100    #Number of samples generated for each matrix
  
  #Condition number of the generated matrix P(W|E,X=1)
  kappas[m] <- kappa(PWEx_t)
  
  U <- matrix(0,n,2)
  W <- matrix(0,n,2)
  X <- matrix(0,n,2)
  Y <- matrix(0,n,2)
  do <- rep(0,N)
  
  PyEx <- rep(0,2)
  PWEx <- matrix(0,2,2)
  
  #Generation of a sample of W in the target domain
  set.seed(1234+n)
  Wnew <- sample(1:2,n,replace=TRUE,prob=QW)
  QW_s <- as.vector(table(factor(Wnew, levels=1:2)))/n
  
  for(j in 1:N){
    set.seed(j)
    for(i in 1:2){
      U[,i] <- sample(1:2,n,replace=TRUE,prob=PUE[,i])
      W[,i] <- sapply(U[,i],substitution)
      X[,i] <- (runif(n)<(0.4+U[,i]/5))*1
      Y[,i] <- (runif(n)<(0.1+X[,i]*0.5+U[,i]/25+W[,i]/50))*1
      
      #Estimation of P(Y=1|E,X=1) and P(W|E,X=1)
      PyEx[i] <- mean(Y[X[,i]==1,i])
      PWEx[,i] <- as.vector(table(factor(W[X[,i]==1,i], levels=1:2)))/sum(X[,i]==1)
    }
    #Estimation of the causal effect
    do[j] <- PyEx%*%solve(PWEx)%*%QW_s
  }
  means[m] <- mean(do)
  confidence[1,m] <- quantile(do,0.025)
  confidence[2,m] <- quantile(do,0.975)
}

#Representation of the estimation of the causal effect depending on the condition number
data <- data.frame(kappas=kappas,means=means,conf1=confidence[1,],conf2=confidence[2,])

ggplot(data=data,aes(x=kappas)) +
  geom_point(aes(y=means)) +
  geom_errorbar(aes(x=kappas,ymin=conf1,ymax=conf2),width=0.5) +
  geom_line(aes(y=effect),color='red',linetype='dashed') +
  theme_bw() +
  labs(x='Condition number of P(W|E,X=1)', y='Q(y|do(x))')


#Relation between condition numbers of P(W|U) and P(W|E,X=1)
K <- 100
kappas1 <- kappas2 <- rep(0,K)

for(k in 1:K){
  set.seed(10+k)
  p <- runif(2)
  PUE <- cbind(c(5,5),c(7,3))/10                           #P(U|E)
  PWU <- cbind(c(p[1],1-p[1]),c(p[2],1-p[2]))/10           #P(W|U) 
  QW <- c(0.4,0.6)                                         #Q(W)
  
  PX1U <-  c(6,8)/10                       #P(X=1|U)
  PX1E <- as.vector(PX1U%*%PUE)            #P(X=1|E)
  PUEX1 <- t(t(PX1U*PUE)/PX1E)             #P(U|E,X=1)
  PWEx_t <- PWU%*%PUEX1                    #P(W|E,X=1)
  
  kappas1[k] <- kappa(PWU)
  kappas2[k] <- kappa(PWEx_t)
}

data_kappa <- data.frame(kappas1, kappas2)

ggplot(data=data_kappa, aes(x=kappas1, y=kappas2)) +
  geom_point() +
  theme_bw() +
  labs(x='Condition number of P(W|U)', y='Condition number of P(W|E,X=1)')

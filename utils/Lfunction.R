#Likelihood function for causal parametrization
Lred <- function(theta_unr,tab_source,tab_target){
  k_Y <- dim(tab_source)[1]
  k_W <- k_U <- dim(tab_source)[2]
  k_E <- dim(tab_source)[3]
  
  #Apply logistic function to obtain parameters in [0,1]
  theta <- 1/(1+exp(-theta_unr))
  
  #P(U|E)
  PUE <- matrix(theta[1:(k_E*(k_U-1))],k_U-1,k_E)
  PUE <- rbind(PUE,1-colSums(PUE))
  
  #P(W|U)
  PWU <- matrix(theta[(k_E*(k_U-1)+1):(k_E*(k_U-1)+k_U*(k_W-1))],k_W-1,k_U)
  PWU <- rbind(PWU,1-colSums(PWU))
  
  #Q(U)
  QU <- theta[(k_E*(k_U-1)+k_U*(k_W-1)+1):((k_E+1)*(k_U-1)+k_U*(k_W-1))]
  QU <- c(QU,1-sum(QU))
  
  #P(doX|U)
  PxU <- theta[((k_E+1)*(k_U-1)+k_U*(k_W-1)+1):((k_E+1)*(k_U-1)+k_U*k_W)]
  
  #P(Y|U,W,x)
  PYUWx <- array(0,c(k_Y,k_U,k_W))
  PYUWx[1:(k_Y-1),,] <- theta[((k_E+1)*(k_U-1)+k_U*k_W+1):((k_E+1)*(k_U-1)+k_U*k_W*k_Y)]
  PYUWx[k_Y,,] <- 1-colSums(PYUWx)
  
  #Q(W)
  QW <- as.numeric(PWU%*%QU)
  
  #Part of the log-likelihood corresponding to the target domain
  L <- -sum(tab_target*log(QW))
  
  PYWxE <- array(NA,c(k_Y,k_W,k_E))
  for(E in 1:k_E){
    #All the values in this loop correspond to a specific domain
    #P(U)
    PU <- PUE[,E]
    
    #P(X=doX)
    Px <- as.numeric(PxU%*%PU)
    
    #P(U|X=doX)
    PUx <- PxU*PU/Px
    
    #P(Y,W,X=doX)
    for(s in 1:k_W){
      PYWxE[,s,E] <- PYUWx[,,s]%*%diag(PUx)%*%t(PWU)[,s]*Px
    }
  }
  
  #Part of the log-likelihood corresponding to the source domain
  L <- L-sum(tab_source*log(PYWxE))

  return(ifelse(is.finite(L),L,1e10))
}

#Function g_y to define the causal estimator
g <- function(pars,k_Y,k_W,k_E,y){
  k_U <- k_W
  
  #P(Y|U,W,X=doX) optimal
  PYUWxo <- array(0,c(k_Y,k_U,k_W))
  PYUWxo[1:(k_Y-1),,] <- pars[((k_E+1)*(k_U-1)+k_U*k_W+1):((k_E+1)*(k_U-1)+k_U*k_W*k_Y)]
  PYUWxo[k_Y,,] <- 1-colSums(PYUWxo)
  
  #P(W|U) optimal
  PWUo <- matrix(pars[(k_E*(k_U-1)+1):(k_E*(k_U-1)+k_U*(k_W-1))],k_W-1,k_U)
  PWUo <- rbind(PWUo,1-colSums(PWUo))
  
  #Q(U) optimal
  QUo <- pars[(k_E*(k_U-1)+k_U*(k_W-1)+1):((k_E+1)*(k_U-1)+k_U*(k_W-1))]
  QUo <- c(QUo,1-sum(QUo))
  
  #Estimation of the causal effect
  estim_effect <- as.numeric(diag(PYUWxo[y,,]%*%PWUo)%*%QUo)
  
  return(estim_effect)
}

#Function h to define the reduced estimator
h <- function(eta,k_W,k_E){
  #Construction of the matrices Q(W), P(W|E,x) and P(y|E,x) from eta
  QW <- c(eta[1:(k_W-1)]/eta[k_W],1-sum(eta[1:(k_W-1)]/eta[k_W]))

  PWEx_num <- matrix(eta[(k_W+1):(k_W+(k_W-1)*k_E)],k_W-1,k_E)
  PWEx_den <- eta[(k_W+k_W*k_E+1):(k_W+(k_W+1)*k_E)]
  PWEx0 <- PWEx_num%*%diag(noNaN(1/PWEx_den))
  PWEx <- rbind(PWEx0,1-colSums(PWEx0))

  PyEx_num <- eta[(k_W+(k_W-1)*k_E+1):(k_W+k_W*k_E)]
  PyEx <- PyEx_num*noNaN(1/PWEx_den)
  
  #Estimation of the causal effect
  if(kappa(PWEx)<1e14){
    estim_effect <- as.numeric(PyEx%*%pseudosolve(PWEx)%*%QW)
  }else{
    estim_effect <- 1
  }
  
  return(estim_effect)
}

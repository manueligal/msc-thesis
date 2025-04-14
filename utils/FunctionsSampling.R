source('misc.R')

#Generation of the matrices that represent the SCM
MechanismGeneration <- function(k_Y,k_X,k_W,k_E,y,doX){
  k_U <- k_W

  PUE <- mat_gen(k_U,k_E)           #P(U|E)
  QU <- mat_gen(k_U,1)              #Q(U)
  PWU <- mat_gen(k_W,k_U)           #P(W|U) 
  PXU <- mat_gen(k_X,k_U)           #P(X|U)
  PYUWX <- array(t(rdirichlet(k_U*k_W*k_X,rep(1,k_Y))),c(k_Y,k_U,k_W,k_X))  #P(Y|U,W,X)
  
  #Other values of interest
  QW <- PWU%*%QU                                #Q(W)
  PxU <- PXU[doX,]                              #P(X=doX|U)
  PxE <- as.vector(PxU%*%PUE)                   #P(X=doX|E)
  PUEx <- t(t(PxU*PUE)/PxE)                     #P(U|E,X=doX)
  PWEx <- PWU%*%PUEx                            #P(W|E,X=doX)
  PyUWx <- PYUWX[y,,,doX]                       #P(Y=y|U,W,X=doX)
  PyUx <- diag(PyUWx%*%PWU)                     #P(Y=y|U,X=doX)
  PyEx <- PyUx%*%PUEx                           #P(Y=y|E,X=doX)
  effect <- PyEx%*%pseudosolve(PWEx)%*%QW       #q(Y=y|do(X=doX))
  
  return(list(PUE=PUE,QU=QU,PWU=PWU,PXU=PXU,PYUWX=PYUWX,PWEx=PWEx,effect=effect))
}

#Generation of a sample of size n from the matrices that describe the mechanism
SampleGeneration <- function(matrixList,n){
  #Extract the matrices from the list
  PUE <- matrixList$PUE
  PWU <- matrixList$PWU
  PXU <- matrixList$PXU
  PYUWX <- matrixList$PYUWX
  QU <- matrixList$QU
  PUE_ext <- cbind(PUE,QU)
  
  #Calculate the dimensions
  k_Y <- dim(PYUWX)[1]
  k_X <- nrow(PXU)
  k_U <- k_W <- nrow(PUE)
  k_E <- ncol(PUE)
  
  #Matrix to store the sample
  Sample <- matrix(NA,n,5)
  colnames(Sample) <- c('E','U','W','X','Y')
  
  #Generate the sample following the SCM
  Sample[,'E'] <- sample(1:(k_E+1),n,replace=TRUE) 
  Sample[,'U'] <- sapply(Sample[,'E'],substitution,P=PUE_ext)
  Sample[,'W'] <- sapply(Sample[,'U'],substitution,P=PWU)
  Sample[,'X'] <- sapply(Sample[,'U'],substitution,P=PXU)
  Sample[,'Y'] <- mapply(substitution2,Sample[,'U'],Sample[,'W'],Sample[,'X'],MoreArgs=list(P=PYUWX))
  
  return(Sample)
}

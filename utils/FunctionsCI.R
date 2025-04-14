library(pracma)

#Calculate the standard deviation for asymptotic normal approximation in the reduced parameterization
ReducedSigma <- function(Sample,y,doX){
  #Calculate the dimensions
  k_Y <- max(Sample[,'Y'])
  k_W <- max(Sample[,'W'])
  k_E <- max(Sample[,'E'])-1
  
  #Matrix with the vectors eta^i
  etas <- matrix(NA,nrow(Sample),k_W+(k_W+1)*k_E)
  
  for(j in 1:(k_W-1)){
    etas[,j] <- ((Sample[,'W']==j)&(Sample[,'E']==(k_E+1)))*1
  }
  
  etas[,k_W] <- (Sample[,'E']==(k_E+1))*1
  
  for(l in 1:k_E){
    for(j in 1:(k_W-1)){
      etas[,k_W+(l-1)*(k_W-1)+j] <- ((Sample[,'W']==j)&(Sample[,'X']==doX)&(Sample[,'E']==l))*1
    }
    etas[,k_W+(k_W-1)*k_E+l] <- ((Sample[,'Y']==y)&(Sample[,'X']==doX)&(Sample[,'E']==l))*1
    etas[,k_W*(k_E+1)+l] <- ((Sample[,'X']==doX)&(Sample[,'E']==l))*1
  }
  
  #Sample covariance matrix divided by n
  Sigma <- cov(etas)/nrow(Sample)

  #Estimator for eta
  eta_hat <- colMeans(etas)
  
  #Gradient of h at eta_hat
  h_grad <- grad(h,eta_hat,k_W=k_W,k_E=k_E)
  
  #Estimate for the standard deviation of the asymptotic normal distribution
  sigma <- as.numeric(sqrt(h_grad%*%Sigma%*%h_grad))

  return(sigma)
}

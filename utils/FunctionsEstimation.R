source('Lfunction.R')

#Estimate the causal effect using the reduced parameterization
ReducedEffectEstimation <- function(Sample,y,doX){
  k_W <- max(Sample[,'W'])
  k_Y <- max(Sample[,'Y'])
  
  #Selection of the observations where X=doX
  SampledoX <- Sample[Sample[,'X']==doX,]
  
  #Estimation of the matrices
  QW <- as.vector(table(factor(Sample[which(Sample[,'E']==(k_E+1)),'W'],levels=1:k_W))/sum(Sample[,'E']==(k_E+1)))
  
  PWEx_n0 <- table(factor(SampledoX[,'W'],levels=1:k_W),SampledoX[,'E'])
  PWEx_n <- PWEx_n0[,-which(colnames(PWEx_n0)==as.character(k_E+1))]
  PWEx <- sweep(PWEx_n,2,colSums(PWEx_n),'/')
  
  PYEx_n0 <- table(factor(SampledoX[,'Y'],levels=1:k_Y),SampledoX[,'E'])
  PYEx_n <- PYEx_n0[,-which(colnames(PYEx_n0)==as.character(k_E+1))]
  PyEx <- sweep(PYEx_n,2,colSums(PYEx_n),'/')[y,]
  
  #Estimation of the causal effect
  if(kappa(PWEx)<1e14){
    estim_effect <- as.numeric(PyEx%*%pseudosolve(PWEx)%*%QW)
  }else{
    estim_effect <- NA
  }
  
  return(list(estim_effect=estim_effect,PyEx=PyEx,PWEx=PWEx,QW=QW))
}

#Estimate the causal effect using the causal parameterization
CausalEffectEstimation <- function(Sample,y,doX,Nseeds){
  #Calculate the dimensions
  k_Y <- max(Sample[,'Y'])
  k_U <- k_W <- max(Sample[,'W'])
  k_E <- max(Sample[,'E'])-1
  
  #Dimension of the parameter that is optimized
  par_dim <- (k_U-1)*(k_E+1)+k_Y*k_U*k_W
  
  #Table of the number of observations n(y,doX,w,e)
  tab_source <- table(Sample[,'Y'],Sample[,'X'],Sample[,'W'],Sample[,'E'])[,doX,,1:k_E]
  tab_target  <- table(Sample[Sample[,'E']==(k_E+1),'W'])
  
  k_max <- max(k_W-1,k_Y-1)
  if(k_max>1){
    lim1 <- -log(2*(k_max-1))
    lim2 <- -log(k_max-1)
  }else{
    lim1 <- 0
    lim2 <-1
  }
  
  theta_unr_opt <- NA
  L_opt <- 1e10
  for(seed in 1:Nseeds){
    opt <- optim(runif(par_dim,lim1,lim2),Lred,method='L-BFGS-B',tab_source=tab_source,tab_target=tab_target,control=list(maxit=5e4))
    if(opt$value<L_opt){
      L_opt <- opt$value
      theta_unr_opt <- opt$par
    }
  }
  
  #Transformation to keep the components in [0,1]
  pars  <- 1/(1+exp(-theta_unr_opt))
  
  #P(U|E) optimal
  PUEo <- matrix(pars[1:(k_E*(k_U-1))],k_U-1,k_E)
  PUEo <- rbind(PUEo,1-colSums(PUEo))
  
  #P(W|U) optimal
  PWUo <- matrix(pars[(k_E*(k_U-1)+1):(k_E*(k_U-1)+k_U*(k_W-1))],k_W-1,k_U)
  PWUo <- rbind(PWUo,1-colSums(PWUo))
  
  #Q(U) optimal
  QUo <- pars[(k_E*(k_U-1)+k_U*(k_W-1)+1):((k_E+1)*(k_U-1)+k_U*(k_W-1))]
  QUo <- c(QUo,1-sum(QUo))
  
  #P(doX|U) optimal
  PxUo <- pars[((k_E+1)*(k_U-1)+k_U*(k_W-1)+1):((k_E+1)*(k_U-1)+k_U*k_W)]
  
  #P(y|U,W,X=doX) optimal
  PYUWxo <- array(0,c(k_Y,k_U,k_W))
  PYUWxo[1:(k_Y-1),,] <- pars[((k_E+1)*(k_U-1)+k_U*k_W+1):((k_E+1)*(k_U-1)+k_U*k_W*k_Y)]
  PYUWxo[k_Y,,] <- 1-colSums(PYUWxo)
  PyUWxo <- PYUWxo[y,,]
  
  #Estimation of the causal effect
  estim_effect <- as.numeric(diag(PyUWxo%*%PWUo)%*%QUo)
  
  return(list(estim_effect=estim_effect,pars=pars,L=L_opt))
}

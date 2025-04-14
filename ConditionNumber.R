library(ggplot2)
library(reshape2)
source('utils/FunctionsSampling.R',chdir=TRUE)
source('utils/FunctionsEstimation.R',chdir=TRUE)

#Generation parameters
n <- 2e4  
M <- 25      
N <- 5    

#Number of categories
k_E <- 2  
k_W <- 2  
k_X <- 2  
k_Y <- 2  

#Study the causal effect q(Y=y|do(X=doX)). The value in () is in binary notation
doX <- (1)+1
y <- (1)+1

data <- data.frame(matrix(NA,M*N,4))
colnames(data) <- c('kappa_theor','kappa_samp','Causal','Reduced')

set.seed(127)

for(m in 1:M){
  matrices <- MechanismGeneration(k_Y,k_X,k_W,k_E,y,doX)
  data[((m-1)*N+1):(N*m),'kappa_theor'] <- kappa(matrices$PWEx)
  
  for(j in 1:N){
    Sample <- SampleGeneration(matrices,n)
    
    #Estimation with the reduced parameterization
    reduced_estim <- ReducedEffectEstimation(Sample,y,doX)
    data[(m-1)*N+j,'Reduced'] <- reduced_estim$estim_effect-matrices$effect
    
    #Condition number of the estimation of P(W|E,X=doX)
    data[N*(m-1)+j,'kappa_samp'] <- kappa(reduced_estim$PWEx)
    
    #Estimation with the causal parameterization
    optimal_effect <- CausalEffectEstimation(Sample,y,doX,3)$estim_effect
    data[(m-1)*N+j,'Causal'] <- optimal_effect-matrices$effect
  }
}

data['sample'] <- 1:nrow(data)
data_melt <- melt(data,id.vars=c('kappa_theor','kappa_samp','sample'),variable.name='method',value.name='error')

#Estimation error as a function of condition number
ggplot(data=data_melt,aes(x=log(kappa_samp),y=abs(error),color=log(kappa_theor))) +
  facet_wrap(~method) +
  geom_point() +
  labs(x=expression('log('~kappa~'('~widehat(P)[n]~'(W|E,x)))'), y='Estimation error',
       color=expression('log('~kappa~'('~'P(W|E,x)))')) +
  theme_bw() +
  theme(legend.direction='vertical',legend.position='right',
        legend.title=element_text(margin=margin(0,0,0,-30,'pt')),
        legend.box.margin=margin(0,0,0,25)) +
  ylim(0,1.5)

#Simulations with one additional domain
k_E <- 3
n <- 2.5e4
data_extra <- data.frame(matrix(NA,M*N,4))
colnames(data_extra) <- c('kappa_theor','kappa_samp','Causal','Reduced')

set.seed(128)

for(m in 1:M){
  matrices <- MechanismGeneration(k_Y,k_X,k_W,k_E,y,doX)
  data_extra[((m-1)*N+1):(N*m),'kappa_theor'] <- kappa(matrices$PWEx)
  
  for(j in 1:N){
    Sample <- SampleGeneration(matrices,n)
    
    #Estimation with the reduced parameterization
    reduced_estim <- ReducedEffectEstimation(Sample,y,doX)
    data_extra[(m-1)*N+j,'Reduced'] <- reduced_estim$estim_effect-matrices$effect
    
    #Condition number of the estimation of P(W|E,X=doX)
    data_extra[N*(m-1)+j,'kappa_samp'] <- kappa(reduced_estim$PWEx)
    
    #Estimation with the causal parameterization
    optimal_effect <- CausalEffectEstimation(Sample,y,doX,3)$estim_effect
    data_extra[(m-1)*N+j,'Causal'] <- optimal_effect-matrices$effect
  }
}

data_extra['sample'] <- 1:nrow(data)
data_extra_melt <- melt(data_extra,id.vars=c('kappa_theor','kappa_samp','sample'),variable.name='method',value.name='error')

#Estimation error as a function of condition number
ggplot(data=data_extra_melt,aes(x=log(kappa_samp),y=abs(error),color=log(kappa_theor))) +
  facet_wrap(~method) +
  geom_point() +
  labs(x=expression('log('~kappa~'('~widehat(P)[n]~'(W|E,x)))'), y='Estimation error',
       color=expression('log('~kappa~'('~'P(W|E,x)))')) +
  theme_bw() +
  theme(legend.direction='vertical',legend.position='right',
        legend.title=element_text(margin=margin(0,0,0,-30,'pt')),
        legend.box.margin=margin(0,0,0,25)) +
  ylim(0,1.5)

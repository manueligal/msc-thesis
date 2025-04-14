library(ggplot2)
source('utils/FunctionsSampling.R',chdir=TRUE)
source('utils/FunctionsEstimation.R',chdir=TRUE)
source('utils/FunctionsCI.R',chdir=TRUE)

#Parameters for generation and bootstrap
sample_sizes <- 10^(3:5)
M <- 20
N <- 5
N_boot <- 100

#Number of categories
k_W <- 2
k_E <- 2
k_Y <- 2
k_X <- 2

#Study the causal effect q(Y=y|do(X=doX)). The value in () is in binary notation
doX <- (1)+1
y <- (1)+1

#Confidence levels
alphas <- c(0.10,0.05,0.01)
method <- c('Bootstrap','Reduced')

#Construction of the dataframe
rows <- length(sample_sizes)*2*length(alphas)
data_b <- data.frame(matrix(NA,rows,5))
colnames(data_b) <- c('n','method','alpha','coverage','length')
data_b['n'] <- rep(sample_sizes,each=2*length(alphas))
data_b['method'] <- rep(method,each=3)
data_b['alpha'] <- as.factor(rep(1-alphas,length(sample_sizes)))

set.seed(123)

for(size in 1:length(sample_sizes)){
  n <- sample_sizes[size]
  
  #Arrays to store if the CI covers the true causal effect and its length for each sample
  covered <- array(NA,c(M*N,3,2))
  length <- array(NA,c(M*N,3,2))
  
  for(m in 1:M){
    matrices <- MechanismGeneration(k_Y,k_X,k_W,k_E,y,doX)
    
    for(j in 1:N){
      Sample <- SampleGeneration(matrices,n)
      
      #Estimation using the reduced estimator
      estim_effect_R <- ReducedEffectEstimation(Sample,y,doX)$estim_effect
      sigma_R <- ReducedSigma(Sample,y,doX)
      
      #Estimation of the standard deviation using bootstrap
      estim_boot <- rep(NA,N_boot)
      for(b in 1:N_boot){
        missing <- 0
        while(missing<(k_E+1)){
          #Generate the bootstrap sample
          Sample_boot <- Sample[sample(1:n,n,replace=TRUE),]
          missing <- length(unique(Sample_boot[Sample_boot[,'X']==doX,'E']))
        }
        #Estimate the causal effect
        estim_boot[b] <- ReducedEffectEstimation(Sample_boot,y,doX)$estim_effect
      }
      #Calculate the standard deviation 
      sigma_boot <- sd(estim_boot)
      
      for(a in 1:length(alphas)){
        alpha <- alphas[a]
        
        #Confidence intervals for bootstrap
        CIsides_boot <- estim_effect_R+c(-1,1)*sigma_boot*qnorm(1-alpha/2)
        covered[(m-1)*N+j,a,1] <- ((matrices$effect>=CIsides_boot[1])&(matrices$effect<=CIsides_boot[2]))*1
        length[(m-1)*N+j,a,1] <- CIsides_boot[2]-CIsides_boot[1]
        
        #Confidence intervals for the reduced parameterization
        CIsides_R <- estim_effect_R+c(-1,1)*sigma_R*qnorm(1-alpha/2)
        covered[(m-1)*N+j,a,2] <- ((matrices$effect>=CIsides_R[1])&(matrices$effect<=CIsides_R[2]))*1
        length[(m-1)*N+j,a,2] <- CIsides_R[2]-CIsides_R[1]
      }
    }
  }
  data_b[(size-1)*2*length(alphas)+(1:3),'coverage'] <- apply(covered[,,1],2,mean,na.rm=TRUE)
  data_b[(size-1)*2*length(alphas)+(4:6),'coverage'] <- apply(covered[,,2],2,mean,na.rm=TRUE)
  data_b[(size-1)*2*length(alphas)+(1:3),'length'] <- apply(length[,,1],2,median,na.rm=TRUE)
  data_b[(size-1)*2*length(alphas)+(4:6),'length'] <- apply(length[,,2],2,median,na.rm=TRUE)
}

ggplot(data_b,aes(x=n,y=coverage)) +
  facet_wrap(~method) +
  geom_point(aes(color=alpha)) +
  geom_line(aes(color=alpha)) +
  theme_bw() +
  theme(panel.spacing = unit(2, "lines")) +
  scale_color_manual(values=c('red','blue','orange')) +
  labs(x='Sample size',y='Coverage',color='Level') +
  scale_x_continuous(breaks=sample_sizes,trans='log10') +
  ylim(0.75,1)

ggplot(data_b,aes(x=n,y=length)) +
  facet_wrap(~method) +
  geom_point(aes(color=alpha)) +
  geom_line(aes(color=alpha)) +
  theme_bw() +
  theme(panel.spacing = unit(2, "lines")) +
  scale_color_manual(values=c('red','blue','orange')) +
  labs(x='Sample size',y='Median CI length',color='Level') +
  scale_x_continuous(breaks=sample_sizes,trans='log10')

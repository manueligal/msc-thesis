library(ggplot2)
library(ggpubr)
source('utils/FunctionsSampling.R',chdir=TRUE)
source('utils/FunctionsEstimation.R',chdir=TRUE)
source('utils/FunctionsCI.R',chdir=TRUE)

#Generation parameters
sample_sizes <- 10^(2:5)
M <- 20
N <- 5

#Number of categories
k_E <- 2  
k_W <- 2  
k_X <- 2  
k_Y <- 2  

#Study the causal effect q(Y=y|do(X=doX)). The value in () is in binary notation
doX <- (1)+1
y <- (1)+1

#Confidence levels
alphas <- c(0.01,0.05,0.10)

#Construction of the dataframe
data_CI <- data.frame(matrix(NA,length(sample_sizes)*length(alphas),4))
colnames(data_CI) <- c('n','alpha','coverage','length')
data_CI['n'] <- rep(sample_sizes,each=length(alphas))
data_CI['alpha'] <- as.factor(rep(1-alphas,length(sample_sizes)))

set.seed(125)

for(size in 1:length(sample_sizes)){
  n <- sample_sizes[size]
  
  #Arrays to store if the CI covers the true causal effect and its length for each sample
  covered <- length <- matrix(NA,M*N,length(alphas))

  for(m in 1:M){
    matrices <- MechanismGeneration(k_Y,k_X,k_W,k_E,y,doX)

    for(j in 1:N){
      Sample <- SampleGeneration(matrices,n)
      
      #Estimation using and reduced estimator
      estim_effect <- ReducedEffectEstimation(Sample,y,doX)$estim_effect
      sigma <- ReducedSigma(Sample,y,doX)
      
      for(a in 1:length(alphas)){
        alpha <- alphas[a]
        CIsides <- estim_effect+c(-1,1)*sigma*qnorm(1-alpha/2)
        covered[(m-1)*N+j,a] <- ((matrices$effect>=CIsides[1])&(matrices$effect<=CIsides[2]))*1
        length[(m-1)*N+j,a] <- CIsides[2]-CIsides[1]
      }
    }
  }
  data_CI[(size-1)*length(alphas)+(1:3),'coverage'] <- apply(covered,2,mean,na.rm=TRUE)
  data_CI[(size-1)*length(alphas)+(1:3),'length'] <- apply(length,2,median,na.rm=TRUE)
}

g1 <- ggplot(data_CI,aes(x=n,y=coverage)) +
  geom_point(aes(color=alpha)) +
  geom_line(aes(color=alpha)) +
  theme_bw() +
  theme(panel.spacing = unit(2, "lines")) +
  scale_color_manual(values=c('red','blue','orange')) +
  labs(x='Sample size',y='Coverage',color='Level') +
  scale_x_continuous(breaks=sample_sizes,trans='log10') +
  ylim(0.75,1)

g2 <- ggplot(data_CI,aes(x=n,y=length)) +
  geom_point(aes(color=alpha)) +
  geom_line(aes(color=alpha)) +
  theme_bw() +
  theme(panel.spacing = unit(2, "lines")) +
  scale_color_manual(values=c('red','blue','orange')) +
  labs(x='Sample size',y='Median CI length',color='Level') +
  scale_x_continuous(breaks=sample_sizes,trans='log10')

ggarrange(g1,g2,nrow=1,legend='right',common.legend=TRUE)

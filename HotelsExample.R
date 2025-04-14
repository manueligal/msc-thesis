library(ggplot2)
source('utils/FunctionsSampling.R',chdir=TRUE)
source('utils/FunctionsEstimation.R',chdir=TRUE)
source('utils/FunctionsCI.R',chdir=TRUE)

#Load the data
train <- read.csv('data/train.csv', header=TRUE)

#Transform the price into a categorical variable
train$price_usd_cat <- as.numeric(cut(train$price_usd,breaks=c(0,75,125,175,225,1e8)))

#Obtain the matrix of the sample of (E,W,X,Y)
sub <- train[,c('prop_id','price_usd_cat','position','click_bool','random_bool')]
sub <- sub[!is.na(sub$price_usd_cat),]
sub$click_bool <- sub$click_bool+1

#Divide into the observational and interventional datasets
obs <- sub[sub$random_bool==0,1:4]
exp <- sub[sub$random_bool==1,1:4]

#Keep only hotels with a minimum sample size
#Observational dataset
hotels_obs <- table(obs$prop_id)
hotels_selection_obs <- as.double(rownames(hotels_obs[hotels_obs>2000]))

#Randomized dataset
hotels_exp <- table(exp$prop_id)
hotels_selection_exp <- as.double(rownames(hotels_exp[hotels_exp>1500]))

#Select the source domains
hotels_kept <- setdiff(hotels_selection_obs,hotels_selection_exp)
obs_selection <- obs[obs$prop_id%in%hotels_kept,]
obs_selection$prop_id <- factor(obs_selection$prop_id)
levels(obs_selection$prop_id) <- 1:length(hotels_kept)
obs_selection$prop_id <- as.numeric(obs_selection$prop_id)

#Selection of parameters
k_E <- length(hotels_kept)
k_Y <- 2
k_X <- 2
k_W <- length(unique(obs$price_usd))
doX <- (1)+1
y <- (1)+1


#Estimate the causal effect for some holdout target domains
results <- data.frame(matrix(NA,length(hotels_selection_exp),12))
colnames(results) <- c('estim','ci1','ci2','exp','ci3','ci4','cond','ci5','ci6')

for(i_id in 1:length(hotels_selection_exp)){
  id <- hotels_selection_exp[i_id]
  
  #Observational dataset corresponding to the target domain
  target <- obs[obs$prop_id==id,]
  target$prop_id <- k_E+1
  Sample <- rbind(obs_selection,target)
  Sample$position <- (Sample$position==1)*1+1
  colnames(Sample) <- c('E','W','X','Y')

  #Estimation of the causal effect
  estim_effect <- ReducedEffectEstimation(Sample,y,doX)$estim_effect
  results[i_id,'estim'] <- estim_effect
  
  #Calculation of the causal effect from the interventional data
  target_exp <- exp[exp$prop_id==id,]
  target_exp$position <- (target_exp$position==1)*1+1
  colnames(target_exp) <- c('E','W','X','Y')
  estim_exp <- mean(target_exp[target_exp$X==doX,'Y']==y)
  n_exp <- length(target_exp[target_exp$X==doX,'Y'])
  results[i_id,'exp'] <- estim_exp

  #Calculation of the causal effect using the conditional distribution
  target$position <- (target$position==1)*1+1
  colnames(target) <- c('E','W','X','Y')
  estim_cond <- mean(target[target$X==doX,'Y']==y)
  n_cond <- length(target[target$X==doX,'Y'])
  results[i_id,'cond'] <- estim_cond

  sigma <- ReducedSigma(Sample,y,doX)
  results[i_id,c('ci1','ci2')] <- estim_effect+c(-1,1)*qnorm(0.975)*sigma
  
  results[i_id,c('ci3','ci4')] <- estim_exp+c(-1,1)*qnorm(0.975)*sqrt(estim_exp*(1-estim_exp)/n_exp)
  
  results[i_id,c('ci5','ci6')] <- estim_cond+c(-1,1)*qnorm(0.975)*sqrt(estim_cond*(1-estim_cond)/n_cond)
}

results_melt1 <- cbind(results[,c(4:6,1:3)],method=rep('Reduced',length(hotels_selection_exp)))
results_melt2 <- cbind(results[,c(4:6,7:9)],method=rep('Conditional',length(hotels_selection_exp)))
colnames(results_melt2)[4:6] <- c('estim','ci1','ci2')
results_melt <- rbind(results_melt1,results_melt2)
results_melt$target <- rep(1:length(hotels_selection_exp),2)
results_melt$method <- factor(results_melt$method,levels=c('Reduced','Conditional'))

ggplot(data=results_melt,aes(x=target)) +
  facet_wrap(~method,nrow=2) +
  geom_point(aes(y=estim,color=method)) +
  geom_errorbar(aes(ymin=(ci1),ymax=(ci2),color=method),width=0.5) +
  geom_point(aes(y=exp,color='Oracle')) +
  geom_errorbar(aes(ymin=(ci3),ymax=(ci4),color='Oracle'),width=0.5) +
  theme_bw() +
  theme(strip.background=element_blank(),strip.text.x=element_blank()) +
  labs(x='Target Domain ID',y=expression(widehat(q)[n]~'(Y=1|do(X=1))')) +
  scale_color_manual('',breaks=c('Reduced','Conditional','Oracle'),values=c('Reduced'='black','Oracle'='red','Conditional'='blue')) +
  coord_cartesian(ylim=c(-0.025,0.425))

#Estimate both effects for different atomic interventions
K <- 20
results_k <- data.frame(matrix(NA,K,6))
colnames(results_k) <- c('estim1','ci1','ci2','exp','ci3','ci4')

#Dataset corresponding to the target domain
id <- hotels_selection_exp[1]
target <- obs[obs$prop_id==id,]
target$prop_id <- k_E+1
target_exp <- exp[exp$prop_id==id,]

for(k in 1:K){
  print(k)
  Sample <- rbind(obs_selection,target)
  Sample$position <- (Sample$position==k)*1+1
  colnames(Sample) <- c('E','W','X','Y')
  
  #Estimation using the reduced parameterization
  estim_effect1 <- ReducedEffectEstimation(Sample,y=2,doX=2)$estim_effect
  results_k[k,'estim1'] <- estim_effect1

  sigma1 <- ReducedSigma(Sample,y,doX)
  results_k[k,c('ci1','ci2')] <- estim_effect1+c(-1,1)*qnorm(0.975)*sigma1
  
  #Calculation of the causal effect from the interventional data
  Sample_exp <- target_exp
  Sample_exp$position <- (Sample_exp$position==k)*1+1
  colnames(Sample_exp) <- c('E','W','X','Y')
  estim_exp <- mean(Sample_exp[Sample_exp$X==doX,'Y']==y)
  n_exp <- length(Sample_exp[Sample_exp$X==doX,'Y'])
  results_k[k,'exp'] <- estim_exp
  
  results_k[k,c('ci3','ci4')] <- estim_exp+c(-1,1)*qnorm(0.975)*sqrt(estim_exp*(1-estim_exp)/n_exp)
}

#We calculate the proportion of clicks across all positions
average_clicks <- mean(target_exp$click_bool==y)

ggplot(data=results_k,aes(x=1:K)) +
  geom_hline(aes(yintercept=average_clicks,linetype='Avg click proportion'),color='red',alpha=0.5) +
  geom_point(aes(y=estim1,color='Reduced')) +
  geom_errorbar(aes(ymin=(ci1),ymax=(ci2),color='Reduced'),width=0.5) +
  geom_point(aes(y=exp,color='Oracle')) +
  geom_errorbar(aes(ymin=(ci3),ymax=(ci4),color='Oracle'),width=0.5) +
  theme_bw() +
  labs(x='Position x',y=expression(widehat(q)[n]~'(Y=1|do(X=x))')) +
  scale_color_manual('',values=c('Reduced'='black','Oracle'='red')) +
  scale_linetype_manual('',values=c('Avg click proportion'='dashed')) +
  coord_cartesian(ylim=c(-0.2,0.5))

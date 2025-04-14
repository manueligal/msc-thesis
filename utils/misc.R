#Script containing auxiliary functions 
library(gtools)

#Create a matrix of size k1xk2 with columns that add up to 1
mat_gen <- function(k1,k2){
  t(rdirichlet(k2,rep(1,k1)))
}

#Obtain values for variables using conditional probability matrices
substitution <- function(i,P){
  k <- dim(P)[1]
  sample(1:k,1,prob=P[,i])
}

substitution2 <- function(i1,i2,i3,P){
  k <- dim(P)[1]
  sample(1:k,1,prob=P[,i1,i2,i3])
}

#Calculate the (right) Moore-Penrose pseudoinverse of a matrix
pseudosolve <- function(A){
  t(A)%*%solve(A%*%t(A))
}

#Transform infinite values into 0
noNaN <- function(x){
  ifelse(is.finite(x),x,0)
}

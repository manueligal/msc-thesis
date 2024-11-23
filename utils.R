#Script containing useful functions 

#Function to generate k numbers that add up to 1 (2 methods)
vec_gen <- function(k,method){
  if(method==1){
    diff(c(0,sort(runif(k-1)),1))
  }else{
    x <- runif(k)
    x/sum(x)
  }
}

#Function to create a matrix of size kxm with columns that add up to 1
mat_gen <- function(k,m){
  replicate(m,vec_gen(k,1))
}

#Function to obtain a sample of W from a sample of U
substitution <- function(i){
  sample(1:k_W,1,prob=PWU[,i])
}

#Calculate the (right) Moore-Penrose pseudoinverse of a matrix
pseudosolve <- function(A){
  t(A)%*%solve(A%*%t(A))
}
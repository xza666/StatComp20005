#' @title The alternating direction method of multipliers (ADMM) algorithm for subgroup analysis with the SCAD penalty
#' @description The ADMM algorithm implements subgroup analysis based on regression models
#' @param x n*p design matrix with p predictors
#' @param y n-dimensional response vector
#' @param lambda the penalty parameter of SCAD
#' @param theta another penalty parameter
#' @param n the sample size
#' @param gamma a parameter that controls the concavity of SCAD
#' @param epsilon the iteration accuracy
#' @return the estimators of intercept, cofficient and the result of subgroup
#' @examples 
#' \dontrun{
#' n <- 100
#' p <- 5
#' x <- matrix(0,n,p)
#' sigma <- matrix(c(rep(0.3,p^2)), p, p)
#' for(i in 1:5) sigma[i,i] <- 1
#' x <- mvrnorm(n,rep(0,p),sigma,empirical = TRUE)
#' eps <- as.matrix(rnorm(n,0,0.25)) 
#' beta0 <- matrix(0,1,p)
#' beta0 <- runif(p,0.5,1)
#' u0 <- as.matrix(2*rbinom(n,1,0.5)-1)   
#' y0 <- matrix(0,n,1)
#' for(i in 1:n) y0[i,] <- u0[i,1]+t(x[i,])%*%beta0+eps[i,]
#' ADMM(x, y0, 0.45, 100, 1, 3, 1e-5)
#' }
#' @export
ADMM <- function(x, y, lambda, n, theta, gamma, epsilon){
  # the smoothly clipped absolute deviation penalty (SCAD) is used in this ADMM algorithm
  
  # define ST function
  ###--------ST--------###
  positive <- function(x){
    if(x >= 0) x <- x
    else x <- 0
    return(x)
  }
  
  ST <- function(t, lambda){
    return(sign(t)*positive(abs(t)-lambda))
  }
  
  ###--------step 1--------###
  # define relevant matrixs
  # these list are used to record iterations in each step
  r <- list()
  mu <- list()
  nu <- list()
  eta <- list()
  beta <- list()
  delte <- list()
  I <- diag(n)
  delta <- matrix(0,n*(n-1)/2,n)
  for(i in 1:(n-1)) {
    for(j in (i+1):n)  { if(i==1) delta[j-i,] <- I[,i]-I[,j] 
    else    delta[(2*n-i)*(i-1)/2+j-i,] <- I[,i]-I[,j] }
  }
  
  beta[[1]] <- solve(t(x)%*%x)%*%t(x)%*%y
  mu[[1]] <- y - x%*%beta[[1]]
  r[[1]] <- c(rep(0,n*(n-1)/2))
  nu[[1]] <- c(rep(0,n*(n-1)/2))
  eta[[1]] <- c(rep(0,n*(n-1)/2))
  delte[[1]] <- c(rep(0,n*(n-1)/2))
  
  k <- 1
  for(i in 1:n){
    for(j in 1:n){
      if(i<j){
        eta[[1]][k] <- mu[[1]][i] - mu[[1]][j]
        k <- k+1
      }
    }
  }
  r[[1]] <- delta%*%mu[[1]] - eta[[1]]
  
  ###--------step 2&3--------###
  m <- 1 # record iteration steps
  Qx <- x%*%solve(t(x)%*%x)%*%t(x)
  # detailed iteration process
  while(m==1||sum(r[[m]]^2) > epsilon){
    # update mu and beta
    mu[[m+1]] <- solve(theta*t(delta)%*%delta+I-Qx)%*%((I-Qx)%*%y+theta*t(delta)%*%(eta[[m]]-theta^(-1)*nu[[m]]))
    beta[[m+1]] <- solve(t(x)%*%x)%*%t(x)%*%(y - mu[[m+1]])
    
    k <- 1
    delte[[m+1]] <- c(rep(0,n*(n-1)/2))
    eta[[m+1]] <- c(rep(0,n*(n-1)/2))
    # update eta
    for(i in 1:n){
      for(j in 1:n){
        if(i < j){
          delte[[m+1]][k] <- mu[[m+1]][i] - mu[[m+1]][j] + theta^(-1)*nu[[m]][k]
          a <- abs(delte[[m+1]][k])
          if(a<=(lambda+lambda/theta)){
            eta[[m+1]][k] <- ST(delte[[m+1]][k],lambda/theta)
          }
          if(a>(lambda+lambda/theta) & a<=(gamma*lambda)){
            eta[[m+1]][k] <- ST(delte[[m+1]][k],gamma*lambda/((gamma-1)*theta))/(1-1/((gamma-1)*theta))
          }
          if(a>(gamma*lambda)){eta[[m+1]][k] <- delte[[m+1]][k]}
          k <- k+1
        }
      }
    }
    
    k <- 1
    # update nu
    nu[[m+1]] <- c(rep(0,n*(n-1)/2))
    for(i in 1:n){
      for(j in 1:n){
        if(i < j){
          nu[[m+1]][k] <- nu[[m]][k]+theta*(mu[[m+1]][i]-mu[[m+1]][j]-eta[[m+1]][k])
          k <- k+1
        }
      }
    }
    r[[m+1]] <- delta%*%mu[[m+1]] - eta[[m+1]]
    m <- m+1
  }
  # the estimators are recorded in a list
  return(list(mu_hat = mu[[m]],beta_hat = beta[[m]],eta_hat=eta[[m]]))
}  

#' @title The mean imputation method for modifying right censoring survival time data
#' @description The mean imputation method for modifying right censoring survival time data with K-M estimator
#' @param y0 the n-dimensional vector for the logarithm of failure time
#' @param n the sample size
#' @param per the censoring rate 
#' @importFrom stats runif
#' @return the modified survival data
#' @examples 
#' \dontrun{
#' n <- 100
#' p <- 5
#' x <- matrix(0,n,p)
#' sigma <- matrix(c(rep(0.3,p^2)), p, p)
#' for(i in 1:5) sigma[i,i] <- 1
#' x <- mvrnorm(n,rep(0,p),sigma,empirical = TRUE)
#' eps <- as.matrix(rnorm(n,0,0.25)) 
#' beta0 <- matrix(0,1,p)
#' beta0 <- runif(p,0.5,1)
#' u0 <- as.matrix(2*rbinom(n,1,0.5)-1)   
#' y0 <- matrix(0,n,1)
#' for(i in 1:n) y0[i,] <- u0[i,1]+t(x[i,])%*%beta0+eps[i,]
#' mean_imputation(y0, 100, 0.1)
#' }
#' @export
mean_imputation <- function(y0, n, per) {
  # the logarithm of the complete survival time data
  # the size of sample
  # the censoring rate
  T <- exp(y0)
  T1 <- T
  a <- 1
  k <- 1:n
  while(a<=n*per)  {   ran1 <- runif(1, 0, sort(T)[0.95*n])
  ran2 <- sample(k,1,replace=TRUE)
  if(ran2==0) next;
  if(ran1<T1[ran2]&&k[ran2]==ran2) { a <- a + 1
  k[ran2] <- 0
  T1[ran2] <- ran1 
  } 
  }
  sdelta <- c(rep(1,n)) # the censoring indicator
  for(i in 1:n) if(T1[i]<T[i]) sdelta[i] <- 0
  
  C <- T1[which(sdelta==0)]  # observed right censoring times
  C1 <- T1[which(sdelta==1)]  # observed true failure times
  if(max(T1)==max(C)) { C1 <- c(C1,max(C))
  C <- C[which(C!=max(C))]
  sdelta[which(T1==max(T1))] <- 1 }  
  # take the biggest observed time as true failure time
  tau <- sort(unique(C))  # different censoring times
  tau1 <- sort(unique(C1)) # different failure times
  
  
  len <- length(tau)  
  Nc <- c(rep(0,len))
  R <- c(rep(0,len)) 
  for(i in 1:len) { R[i] <- length(which(C==tau[i]))
  Nc[i] <- length(which(T1>=tau[i]))-length(which(C1==tau[i])) }
  
  len1 <- length(tau1)
  Nc1 <- c(rep(0,len1))
  R1 <- c(rep(0,len1))
  for(i in 1:len1) { R1[i] <- length(which(C1==tau1[i]))
  Nc1[i] <- length(which(T1>=tau1[i])) }
  
  # define Sc function
  # the K-M estimator of the survival function of the right censoring variable
  Sc <- function(x) { g <- which(tau<=x)
  if(length(g)==0) return(1)
  else return(prod(1-R[1:max(g)]/Nc[1:max(g)])) }
  
  # define S function
  # the K-M estimator of the survival function of the true failure time variable
  S <- function(x) { g <- which(tau1<=x)
  if(length(g)==0) return(1)
  else return(prod(1-R1[1:max(g)]/Nc1[1:max(g)])) }
  
  tau2 <- c(0,tau1)
  triS <- c(rep(0,len1))
  for(i in 1:len1) triS[i] <- S(tau2[i])-S(tau2[i+1]) 
  
  # modify the observed survival times by the mean imputation method
  y <- matrix(0,n,1)
  for(i in 1:n) { if(sdelta[i]==1) y[i] <- y0[i]
  else y[i] <- sum(log(tau1[which(tau1>T1[i])])*triS[which(tau1>T1[i])])/S(T1[i]) }
  return(list(y_hat=y))
}

#' @title A method for solving an interesting parking problem
#' @description This method can solve the parking problem in fixed interval 
#' @param x the length of parking interval
#' @param n numbers of independent repeated trials 
#' @importFrom stats runif
#' @return numbers of parking cars and the average of parking cars in unit length
#' @examples 
#' \dontrun{
#' parking_solution(1000, 1000)
#' }
#' @export
parking_solution <- function(x, n) {
  # x: the parking interval
  # the number of simulations 
  p1 <- function(a=0, b=x, c=0){
    while(TRUE){
      if(b-a<1) return(c)
      t<- runif(1,a,b)    # give the location of the first car
      if(t+1>b) next
      else {c <- c+1 ; break}
    }
    a1 <- a
    a2 <- t+1
    b1 <- t
    b2 <- b
    c <- p1(a1,b1,c)
    c <- p1(a2,b2,c)
    return(c)
  }
  d <- numeric()
  i <- 0
  while(i<=n){
    d <- c(d,p1(0,x,0))
    i <- i+1
  }
  m1 <- mean(d)
  lim1 <- m1/x
  return(list(M_x=m1, lim_M_x=lim1))
}


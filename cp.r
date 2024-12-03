###Coverage probability for the DMC, LSUR,HSUR,LHSUR
#######################################
library(MASS);library(statmod);library(coda);library(tidyverse)
library(bayesm);library(MCMCpack);library(ggplot2)
###########
# Data generation
# Generates samples from the multivariate normal distribution
prnorm <- function(mu,A,n,p) {
  U <- svd(A)$u
  V <- svd(A)$v
  D <- diag(sqrt(svd(A)$d))
  B <- U %*% D %*% t(V)
  w <- matrix(0,n,p)
  for (i in 1:n){w[i,] <- mu + B%*%rnorm(p)}
  return(w)
}

NN=100#  repeat 100 times
Dcov11=Dcov12=Dcov13=Dcov14=Dcov15=Dcov16=Dcov17=Dcov18=Dcov19=Dcov110=Dcov111=Dcov112=Dcov113=Dcov114=Dcov115=0
Dcov21=Dcov22=Dcov23=Dcov24=Dcov25=Dcov26=Dcov27=Dcov28=Dcov29=Dcov210=Dcov211=Dcov212=Dcov213=Dcov214=Dcov215=Dcov216=0

Lcov11=Lcov12=Lcov13=Lcov14=Lcov15=Lcov16=Lcov17=Lcov18=Lcov19=Lcov110=Lcov111=Lcov112=Lcov113=Lcov114=Lcov115=0
Lcov21=Lcov22=Lcov23=Lcov24=Lcov25=Lcov26=Lcov27=Lcov28=Lcov29=Lcov210=Lcov211=Lcov212=Lcov213=Lcov214=Lcov215=Lcov216=0

Hcov11=Hcov12=Hcov13=Hcov14=Hcov15=Hcov16=Hcov17=Hcov18=Hcov19=Hcov110=Hcov111=Hcov112=Hcov113=Hcov114=Hcov115=0
Hcov21=Hcov22=Hcov23=Hcov24=Hcov25=Hcov26=Hcov27=Hcov28=Hcov29=Hcov210=Hcov211=Hcov212=Hcov213=Hcov214=Hcov215=Hcov216=0

LHcov11=LHcov12=LHcov13=LHcov14=LHcov15=LHcov16=LHcov17=LHcov18=LHcov19=LHcov110=LHcov111=LHcov112=LHcov113=LHcov114=LHcov115=0
LHcov21=LHcov22=LHcov23=LHcov24=LHcov25=LHcov26=LHcov27=LHcov28=LHcov29=LHcov210=LHcov211=LHcov212=LHcov213=LHcov214=LHcov215=LHcov216=0
for(j in 1:NN){
  n<-200#sample size
  p<-15#10
  ###标准化设计矩阵X_1 
  design_matrix1 <-matrix(runif(n*p,-4,4),ncol=p)   
  # 计算每列的均值
  column_means1 <- colMeans(design_matrix1)
  # 中心化设计矩阵
  centered_matrix1 <- sweep(design_matrix1, 2, column_means1)
  # 计算协方差矩阵
  cov_matrix1 <- cov(centered_matrix1)
  cov_matrix1
  # 特征值分解
  eigen_decomp1 <- eigen(cov_matrix1)
  eigen_values1 <- eigen_decomp1$values
  eigen_vectors1 <- eigen_decomp1$vectors
  # 标准正交化设计矩阵
  orthogonal_matrix1 <- centered_matrix1 %*% eigen_vectors1
  # 归一化
  X_1<-sweep(orthogonal_matrix1, 2, sqrt(eigen_values1), FUN="/")
  
  ##标准化设计矩阵X_2
  design_matrix2 <- matrix(runif(n*p,0,10),ncol=p) 
  # 计算每列的均值
  column_means2 <- colMeans(design_matrix2)
  # 中心化设计矩阵
  centered_matrix2 <- sweep(design_matrix2, 2, column_means1)
  # 计算协方差矩阵
  cov_matrix2 <- cov(centered_matrix2)
  cov_matrix2
  # 特征值分解
  eigen_decomp2 <- eigen(cov_matrix2)
  eigen_values2 <- eigen_decomp2$values
  eigen_vectors2 <- eigen_decomp2$vectors
  # 标准正交化设计矩阵
  orthogonal_matrix2 <- centered_matrix2 %*% eigen_vectors1
  # 归一化
  X_2 <- sweep(orthogonal_matrix2, 2, sqrt(eigen_values2), FUN="/")
  B1 <- c(c(3,3,3),rep(0, p-3)) 
  B2 <- c(c(-3,-3,-3),rep(0, p-3))
  #B1 <- c(c(3,5,7),rep(0, p-3)) 
  #B2 <- c(c(2,4,6),rep(0, p-3))
  SS <-matrix(c(1,0,0,1),2,2)
  er <- prnorm(rep(0,len=2),SS,n,2)
  Y_1 <- X_1%*%B1+er[,1]
  Y_2 <- X_2%*%B2+er[,2]
  rho=SS[1,2]/sqrt(SS[1,1]*SS[2,2])
  Niter=10000
  
  ###DMC方法：
  # n=200
  X <- matrix(0,ncol=length(B1)+length(B2),nrow=2*n)
  X[1:n,1:length(B1)] <- X_1
  X[(n+1):(2*n),(length(B1)+1):(length(B1)+length(B2))]<-X_2
  y <- c(Y_1,Y_2)
  
  hat_B1 <- solve(t(X_1)%*%X_1)%*%t(X_1)%*%Y_1
  n1 <- n-p
  n2 <- n-p-1
  m=6000
  omiga=matrix(0,m,4)
  #beta_rho=matrix(0,m,7)
  DMC_b1=matrix(0,m,length(B1))
  DMC_b2=matrix(0,m,length(B2)+1)
  ####开始计算DMC模拟时间
  timestart<-Sys.time();
  for(i in 1:m){
    S <- matrix(0,2,2)
    S[1,1] <- rinvgamma(1,shape=n1/2,scale=sum((Y_1-X_1%*%hat_B1)^2)/2)
    Beta1 <- mvrnorm(n=1,hat_B1,S[1,1]*solve(t(X_1)%*%X_1))
    Z2    <- cbind(X_2,Y_1-X_1%*%Beta1)
    hat_B2 <- solve(t(Z2)%*%Z2)%*%t(Z2)%*%Y_2
    S[2,2]<-rinvgamma(1,shape=n2/2,scale=sum((Y_2-Z2%*%hat_B2)^2)/2)
    Beta2 <- mvrnorm(n=1,hat_B2,S[2,2]*solve(t(Z2)%*%Z2))
    S[1,2] <- S[2,1] <- Beta2[3]*S[1,1]
    S[2,2] <- S[2,2]+(Beta2[3])^2*S[1,1]
    O <- kronecker(solve(S),diag(1,n))
    B <- solve(t(X)%*%O%*%X)%*%t(X)%*%O%*%y
    B <- mvrnorm(n=1,B,solve(t(X)%*%O%*%X))
    omiga[i,]=c(S)
    DMC_b1[i,]=Beta1
    DMC_b2[i,]=Beta2
  }
  timeend<-Sys.time()
  runningtime<-timeend-timestart
  print(runningtime)###############DMC时间
  
  
  ##方法1 Lasso方法
  updateBeta1<-function(Y_1,X_1,Sig_1_2,
                        Lambda_1)
  {  La <-diag(Lambda_1)
  B=t(X_1)%*%X_1+La
  bhat=solve(t(X_1)%*%X_1)%*%t(X_1)%*%Y_1
  m.post <- solve(B)%*%(t(X_1)%*%X_1)%*%bhat
  sig.post <- solve(B)*Sig_1_2
  return(mvrnorm(1,m.post,sig.post))
  }
  updateBeta2<-function(Y_2,X_2,Sig_2_2,
                        Lambda_2, Y_1, X_1, Beta1)
  {
    Z1=Y_1-X_1%*%Beta1
    Z_2=cbind(X_2,Z1 )
    La2 <-diag(Lambda_2)
    B=t(Z_2)%*%Z_2+La2
    b2hat=solve(t(Z_2)%*%Z_2)%*%t(Z_2)%*%Y_2
    m.post <-solve(B)%*%(t(Z_2)%*%Z_2)%*%b2hat
    sig.post=solve(B)*Sig_2_2
    return(mvrnorm(1,m.post,sig.post))
  }
  updateSig_1_2<-function(Y_1,X_1,
                          Lambda_1, Beta1)
  {
    La1=diag(Lambda_1)
    n_0= length(Y_1)+length(Beta1)+1-1
    S1=t(Y_1-X_1%*%Beta1)%*%(Y_1-X_1%*%Beta1)/2
    S2=t(Beta1)%*%solve(La1)%*%Beta1/2
    a.post <- n_0/2
    b.post <- S1+S2
    return(1/rgamma(1,a.post,b.post))
  }
  
  updateSig_2_2 <- function(Y_2,X_2,Lambda_2, Beta2, Y_1,X_1,Beta1){
    La2 <-diag(Lambda_2)
    n_0= length(Y_2)+length(Beta2)+2-1
    Z_2=cbind(X_2, Y_1-X_1%*%Beta1)
    S1=t(Y_2-Z_2%*%Beta2)%*%(Y_2-Z_2%*%Beta2)/2
    S2=t(Beta2)%*%solve(La2)%*%Beta2/2
    a.post <- n_0/2
    b.post <- S1+S2
    return(1/rgamma(1,a.post,b.post))
  }
  updatelambda01 <- function(a1,s1,Beta1,Lambda_1)
  {  s <-sum(Lambda_1)
  a.post=a1*(length(Beta1)+1-1)-1#修改2024.6.17
  b.post=s1*(length(Beta1)+1-1)+s/2
  return(rgamma(1,a.post,b.post))
  }
  updatelambda02 <- function(a2,s2,
                             Beta2,Lambda_2)
  {  s <-sum(Lambda_2)
  a.post=a2*(length(Beta2)+2-1)-1
  b.post=s2*(length(Beta2)+2-1)+s/2
  return(rgamma(1,a.post,b.post))
  }
  #####逆高斯了!!!!
  updatelambda1<-function(Beta1,lambda01,Sig_1_2)
  {  
    a.post=sqrt(lambda01*Sig_1_2)/sqrt(Beta1*Beta1)
    b.post=lambda01/2
    return(rinvgauss(1,a.post,b.post))#高斯#修改2024.6.17
  }
  updatelambda2<-function(Beta2,lambda02,Sig_2_2)
  {  
    a.post=sqrt(lambda02*Sig_2_2)/sqrt(Beta2*Beta2)
    b.post=lambda02/2
    return(rinvgauss(1,a.post,b.post))#高斯#修改2024.6.17
  }
  #Number of iterations
  Niter=10000
  #Creating matrices that will store MCMC samples
  Beta1.out <- matrix(0, Niter, length(B1) )
  Beta2.out <- matrix(0, Niter,length(B2)+2-1)
  Sig1.out <- matrix(0,  Niter)
  Sig2.out <- matrix(0, Niter)
  Lambda01.out <- matrix(0, Niter )
  Lambda02.out <- matrix(0, Niter)
  ################################################
  Lambda1.out <- matrix(0, Niter, length(B1) )
  Lambda2.out <- matrix(0, Niter, length(B2)+2-1)
  ###################################
  #Initial Values
  Beta1.out[1,]<-rep(1,length(B1))
  Beta2.out[1,]<-rep(1,length(B2)+2-1)
  Sig1.out[1]<-1
  Sig2.out[1]<-1
  #Tau.out[1]<-1
  Lambda01.out[1]<-1
  Lambda02.out[1]<-1
  #######################################
  Lambda1.out[1,]<-rep(1,length(B1))
  Lambda2.out[1,]<-rep(1,length(B2)+2-1)
  a1=a2=1
  s1=s2=1.78
  #Gibbs sampler for LSUR
  for(i in 2:Niter){
    Beta1.out[i,] <- updateBeta1(Y_1,X_1,Sig1.out[i-1],
                                 Lambda1.out[i-1,])
    Beta2.out[i,]<-updateBeta2(Y_2,X_2,Sig2.out[i-1],
                               Lambda2.out[i-1,],Y_1, X_1,Beta1.out[i,] )
    Sig1.out[i]<- updateSig_1_2(Y_1,X_1,
                                Lambda1.out[i-1,], Beta1.out[i,])
    Sig2.out[i]<- updateSig_2_2(Y_2,X_2, Lambda2.out[i-1,],
                                Beta2.out[i,], Y_1,X_1,Beta1.out[i,])
    Lambda01.out[i]<-updatelambda01(a1,s1,Beta1.out[i,],
                                    Lambda1.out[i-1,])
    Lambda02.out[i]<-updatelambda02(a2,s2,Beta2.out[i,],
                                    Lambda2.out[i-1,])
    Lambda1.out[i,]<-updatelambda1(Beta1.out[i,],
                                   Lambda01.out[i],Sig1.out[i])
    Lambda2.out[i,]<-updatelambda2(Beta2.out[i,],
                                   Lambda02.out[i],
                                   Sig2.out[i])
  }
  #Burn-in phase
  Nburn <- 4000
  
  b1 <- Beta1.out[-(1:Nburn),]
  b2 <- Beta2.out[-(1:Nburn),]
  ##方法1 结束
  ###############################
  ##方法2：HSUR
  updateBeta1 <- function(Y_1,X_1,Sig_1_2,
                          Tau_1, kapp_1)
  {
    La <-diag(kapp_1)
    A=t(X_1)%*%X_1+solve(Tau_1*La)
    bhat=solve(t(X_1)%*%X_1)%*%t(X_1)%*%Y_1
    sig.post <- solve(A)*Sig_1_2
    m.post <- solve(A)%*% (t(X_1)%*%X_1)%*%bhat
    return(mvrnorm(1,m.post,sig.post))
  }
  updateBeta2 <- function(Y_2,X_2,Sig_2_2,Tau_2,
                          kapp_2, Y_1, X_1, Beta1)
  {
    Z1=Y_1-X_1%*%Beta1
    Z_2=cbind(X_2,Z1)
    La2 <-diag(kapp_2)
    A_2=t(Z_2)%*%Z_2+solve(Tau_2*La2)
    b2hat=solve(t(Z_2)%*%Z_2)%*%t(Z_2)%*%Y_2
    sig.post <- solve(A_2)*Sig_2_2
    m.post <- solve(A_2)%*%(t(Z_2)%*%Z_2)%*%b2hat
    return(mvrnorm(1,m.post,sig.post))
  }
  updateSig_1_2 <- function(Y_1,X_1,Tau_1,
                            kapp_1, Beta1)
  {
    La1=diag(kapp_1)
    n_0= length(Y_1)+length(Beta1)+1-1#改3.30
    S1=t(Y_1-X_1%*%Beta1)%*%(Y_1-X_1%*%Beta1)/2
    S2=t(Beta1)%*%solve(Tau_1*La1)%*%Beta1/2
    a.post <- n_0/2
    b.post <- S1+S2
    return(1/rgamma(1,a.post,b.post))
  }
  updateSig_2_2 <- function(Y_2,X_2,Tau_2, kapp_2,
                            Beta2, Y_1,X_1,Beta1)
  {
    La2 <-diag(kapp_2)
    n_0= length(Y_2)+length(Beta2)+2-1#改3.20
    Z_2=cbind(X_2, Y_1-X_1%*%Beta1)
    S1=t(Y_2-Z_2%*%Beta2)%*%(Y_2-Z_2%*%Beta2)/2
    S2=t(Beta2)%*%solve(Tau_2*La2)%*%Beta2/2
    a.post <- n_0/2
    b.post <- S1+S2
    return(1/rgamma(1,a.post,b.post))
  }
  updatekapp1<-function(Beta1,Tau_1,Sig_1_2,Nu_1)
  {
    b.post=Beta1*Beta1/(2*Tau_1*Sig_1_2)+1/Nu_1
    return(1/rgamma(1,1,b.post))
  }
  updatekapp2<-function(Beta2,Tau_2,Sig_2_2,Nu_2)
  {
    b.post=Beta2*Beta2/(2*Tau_2*Sig_2_2)+1/Nu_2
    return(1/rgamma(1,1,b.post))
  }
  updateNu1<-function(kapp_1)
  {
    Nu_1=1/rgamma(1,1,1/kapp_1+1)
    return(Nu_1)
  }
  updateNu2<-function(kapp_2)
  {
    Nu_2=1/rgamma(1,1,1/kapp_2+1)
    return(Nu_2)
  }
  updatexi1<-function(Tau_1,Beta1)
  {   a=2*(length(Beta1)+1-1)-1
  b=(length(Beta1)+1-1)*(1/Tau_1+1)
  xi=1/rgamma(1,a,b)
  return(xi)
  }
  updatexi2<-function(Tau_2,Beta2)
  {  a=2*(length(Beta2)+2-1)-1
  b=(length(Beta2)+2-1)*(1/Tau_2+1)
  xi2=1/rgamma(1,a,b)
  return(xi2)
  }
  updateTau1 <- function(Beta1,kapp_1,
                         Sig_1_2, xi_1)
  {
    La1 <-diag(kapp_1)
    n_0=2*(length(Beta1)+1-1)-1
    S1=t(Beta1)%*%solve(Sig_1_2*La1)%*%Beta1/2+(length(Beta1)+1-1)/xi_1
    a.post=n_0
    b.post=S1 
    return(1/rgamma(1,a.post,b.post))
  }
  updateTau2 <- function(Beta2,kapp_2,
                         Sig_2_2, xi_2)
  {
    La2 <-diag(kapp_2)
    n_0=2*(length(Beta2)+2-1)-1
    S1=t(Beta2)%*%solve(Sig_2_2*La2)%*%Beta2/2+(length(Beta2)+2-1)/xi_2
    a.post=n_0
    b.post=S1 
    return(1/rgamma(1,a.post,b.post))
  }
  #Creating matrices that will store MCMC samples
  Beta1.out <- matrix(0, Niter, length(B1) )
  Beta2.out <- matrix(0, Niter,length(B2)+2-1)
  
  Sig1.out <- matrix(0,  Niter)
  Sig2.out <- matrix(0, Niter)
  
  Tau1.out <- matrix(0,  Niter)#
  Tau2.out <- matrix(0,  Niter)#
  
  kapp1.out <- matrix(0, Niter, length(B1) )
  kapp2.out <- matrix(0, Niter, length(B2)+2-1)
  
  nu1.out <- matrix(0, Niter, length(B1))
  nu2.out <- matrix(0, Niter, length(B2)+2-1)
  
  xi1.out <- matrix(0, Niter)
  xi2.out <- matrix(0, Niter)
  
  #Initial Values
  Beta1.out[1,]<-rep(1,length(B1))
  Beta2.out[1,]<-rep(1,length(B2)+2-1)
  Sig1.out[1]<-1
  Sig2.out[1]<-1
  Tau1.out[1]<-1
  Tau2.out[1]<-1
  
  kapp1.out[1,]<-rep(1,length(B1))
  kapp2.out[1,]<-rep(1,length(B2)+2-1)
  
  nu1.out[1,]<-rep(1,length(B1))
  nu2.out[1,]<-rep(1,length(B2)+2-1)
  
  xi1.out[1]<-1
  xi2.out[1]<-1
  #########################################3
  #Gibbs sampler
  for(i in 2:Niter){
    Beta1.out[i,] <- updateBeta1(Y_1,X_1,Sig1.out[i-1],
                                 Tau1.out[i-1],kapp1.out[i-1,])
    Beta2.out[i,]<-updateBeta2(Y_2,X_2,Sig2.out[i-1],
                               Tau2.out[i-1],kapp2.out[i-1,],Y_1, X_1,Beta1.out[i,] )
    Sig1.out[i]<- updateSig_1_2(Y_1,X_1,Tau1.out[i-1],
                                kapp1.out[i-1,], Beta1.out[i,])
    Sig2.out[i]<- updateSig_2_2(Y_2,X_2,Tau2.out[i-1], 
                                kapp2.out[i-1,], Beta2.out[i,],
                                Y_1,X_1,Beta1.out[i,])
    kapp1.out[i,]<-updatekapp1(Beta1.out[i,],
                               Tau1.out[i-1],Sig1.out[i],nu1.out[i-1,])
    kapp2.out[i,]<-updatekapp2(Beta2.out[i,],Tau2.out[i-1],
                               Sig2.out[i],nu2.out[i-1,])
    nu1.out[i,]<-updateNu1(kapp1.out[i,])
    nu2.out[i,]<-updateNu2(kapp2.out[i,])
    xi1.out[i]<-updatexi1(Tau1.out[i-1],Beta1.out[i,])
    xi2.out[i]<-updatexi2(Tau2.out[i-1],Beta2.out[i,])
    
    Tau1.out[i]<-updateTau1(Beta1.out[i,],
                            kapp1.out[i,],Sig1.out[i],xi1.out[i])
    Tau2.out[i]<-updateTau2(Beta2.out[i,],
                            kapp2.out[i,],Sig2.out[i],xi2.out[i])
  }
  ################################
  #Burn-in phase
  Nburn <- 4000
  Bb1 <- Beta1.out[-(1:Nburn),]
  Bb2 <- Beta2.out[-(1:Nburn),]    
  
  ###方法3 LHSUR
  updateBeta1<-function(Y_1,X_1,Sig_1_2,
                        Lambda_1)
  {  La <-diag(Lambda_1)
  B=t(X_1)%*%X_1+solve(La)#6.1888888
  bhat=solve(t(X_1)%*%X_1)%*%t(X_1)%*%Y_1
  m.post <- solve(B)%*%(t(X_1)%*%X_1)%*%bhat
  sig.post <- solve(B)*Sig_1_2
  return(mvrnorm(1,m.post,sig.post))
  }
  updateBeta2 <- function(Y_2,X_2,Sig_2_2,
                          Lambda_2, Y_1, X_1, Beta1,Tau_2, Kapp_2)
  {
    Z1=Y_1-X_1%*%Beta1
    Z_2=cbind(X_2,Z1 )
    #La2=diag( cbind(Lambda_2,Tau_2*Kapp_2))
    La2=diag(c(c(Lambda_2),Tau_2*Kapp_2))
    B=t(Z_2)%*%Z_2+solve(La2)
    b2hat=solve(t(Z_2)%*%Z_2)%*%t(Z_2)%*%Y_2
    m.post <-solve(B)%*%(t(Z_2)%*%Z_2)%*%b2hat
    sig.post=solve(B)*Sig_2_2
    return(mvrnorm(1,m.post,sig.post))
  }
  updateSig_1_2 <- function(Y_1,X_1,
                            Lambda_1, Beta1)
  {
    La1=diag(Lambda_1)
    n_0= length(Y_1)+length(Beta1)+1-1
    S1=t(Y_1-X_1%*%Beta1)%*%(Y_1-X_1%*%Beta1)/2
    S2=t(Beta1)%*%solve(La1)%*%Beta1/2
    a.post <- n_0/2
    b.post <- S1+S2
    return(1/rgamma(1,a.post,b.post))
  }
  
  updateSig_2_2 <- function(Y_2,X_2,
                            Lambda_2, Beta2,
                            Y_1,X_1,Beta1,Tau_2, Kapp_2)
  {
    La2=diag(c(c(Lambda_2),Tau_2*Kapp_2))
    n_0= length(Y_2)+length(Beta2)+2-1
    Z_2=cbind(X_2, Y_1-X_1%*%Beta1)
    S1=t(Y_2-Z_2%*%Beta2)%*%(Y_2-Z_2%*%Beta2)/2
    S2=t(Beta2)%*%solve(La2)%*%Beta2/2
    a.post <- n_0/2
    b.post <- S1+S2
    return(1/rgamma(1,a.post,b.post))
  }
  updatelambda01 <- function(a1,s1,Beta1,Lambda_1)
  {  s <-sum(Lambda_1)
  a.post=a1+length(Beta1)
  b.post=s1+s/2
  return(rgamma(1,a.post,b.post))
  }
  updatelambda02 <- function(a2,s2,
                             Beta2,Lambda_2)
  {  s <-sum(Lambda_2)
  a.post=a2+length(Beta2)
  b.post=s2+s/2
  return(rgamma(1,a.post,b.post))
  }
  #####逆高斯了!!!!
  updatelambda1<-function(Beta1,lambda01,Sig_1_2)
  {
    a.post=sqrt(lambda01*Sig_1_2)/sqrt(Beta1*Beta1)
    b.post=lambda01/2
    return(1/rinvgauss(1,a.post,b.post))####高斯
  }
  
  updatelambda2<-function(Beta2,lambda02,Sig_2_2)
  { 
    a.post=sqrt(lambda02*Sig_2_2)/sqrt(Beta2*Beta2)
    b.post=lambda02/2
    return(1/rinvgauss(1,a.post,b.post))####高斯
  }
  ###################
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  updatekapp2<-function(Beta2,Tau_2,Sig_2_2,Nu_2)
  {
    b.post=(Beta2*Beta2)/(2*Tau_2*Sig_2_2)+1/Nu_2
    return(1/rgamma(1,1,b.post))
  }
  updateNu2<-function(Kapp_2)
  {
    Nu_2=1/rgamma(1,1,1/Kapp_2+1)
    return(Nu_2)
  }
  updatexi2<-function(Tau_2)
  {
    xi2=1/rgamma(1,1,(1/Tau_2+1))
    return(xi2)  }
  updateTau2 <- function(Beta2,Kapp_2,
                         Sig_2_2, xi_2)
  {  
    S1=(Beta2*Beta2)/(Sig_2_2%*%Kapp_2%*%2)+1/xi_2
    a.post=1
    b.post=S1 
    return(1/rgamma(1,a.post,b.post))
  }
  ################
  #Creating matrices that will store MCMC samples
  Beta1.out <- matrix(0, Niter, length(B1) )
  Beta2.out <- matrix(0, Niter,length(B2)+1)
  
  Sig1.out <- matrix(0,  Niter)
  Sig2.out <- matrix(0, Niter)
  
  Lambda01.out <- matrix(0, Niter )
  Lambda02.out <- matrix(0, Niter)
  
  Lambda1.out <- matrix(0, Niter, length(B1))
  Lambda2.out <- matrix(0, Niter, length(B2))
  
  Tau2.out <- matrix(0,  Niter)
  kapp2.out <- matrix(0, Niter)
  
  nu2.out <- matrix(0, Niter)
  xi2.out <- matrix(0, Niter)
  ###################################
  #Initial Values
  Beta1.out[1,]<-rep(1,length(B1))
  Beta2.out[1,]<-rep(1,length(B2)+1)
  Sig1.out[1]<-1
  Sig2.out[1]<-1
  
  Lambda01.out[1]<-1
  Lambda02.out[1]<-1
  #######################################
  Lambda1.out[1,]<-rep(1,length(B1))
  Lambda2.out[1,]<-rep(1,length(B2))
  Tau2.out[1]<-1
  kapp2.out[1]<-1
  nu2.out[1]<-1
  xi2.out[1]<-1
  a1=a2=1
  s1=s2=1.78
  ##方法3  Gibbs sampler for LHSUR
  for(i in 2:Niter){
    Beta1.out[i,] <- updateBeta1(Y_1,X_1,Sig1.out[i-1],
                                 Lambda1.out[i-1,])
    Beta2.out[i,]<-updateBeta2(Y_2,X_2,Sig2.out[i-1],
                               Lambda2.out[i-1,],
                               Y_1, X_1,
                               Beta1.out[i,],
                               Tau2.out[i-1], kapp2.out[i-1])
    
    Sig1.out[i]<- updateSig_1_2(Y_1,X_1, Lambda1.out[i-1,],
                                Beta1.out[i,])
    Sig2.out[i]<- updateSig_2_2(Y_2,X_2,Lambda2.out[i-1,],
                                Beta2.out[i,],
                                Y_1,X_1,Beta1.out[i,],
                                Tau2.out[i-1],kapp2.out[i-1])
    
    Lambda01.out[i]= updatelambda01(a1,s1,
                                    Beta1.out[i,],
                                    Lambda1.out[i-1,])
    Lambda02.out[i]= updatelambda01(a2,s2,
                                    Beta2.out[i,],
                                    Lambda2.out[i-1,])
    
    Lambda1.out[i,]= updatelambda1(Beta1.out[i,],Lambda01.out[i],
                                   Sig1.out[i])
    Lambda2.out[i,]= updatelambda1(Beta2.out[i,],Lambda02.out[i],
                                   Sig2.out[i])
    ################马蹄
    kapp2.out[i]<-updatekapp2(Beta2.out[i,length(B2)+1],
                              Tau2.out[i-1],
                              Sig2.out[i],nu2.out[i-1])
    nu2.out[i]<-updateNu2(kapp2.out[i])
    xi2.out[i]<-updatexi2(Tau2.out[i-1])
    
    Tau2.out[i]<-updateTau2(Beta2.out[i,length(B2)+1],
                            kapp2.out[i],Sig2.out[i],xi2.out[i])
  }
  #Burn-in phase
  Nburn <- 4000
  ########################################
  
  
  Cb1 <- Beta1.out[-(1:Nburn),]
  Cb2 <- Beta2.out[-(1:Nburn),]
  
  
  ################结果0: DMC
  DCI1=apply(DMC_b1,2,quantile, probs=c(0.025,0.975) )
  DCI2=apply(DMC_b2,2,quantile, probs=c(0.025,0.975) )
  
  ###############结果1
  LCI1=apply(b1,2,quantile, probs=c(0.025,0.975) )
  LCI2=apply(b2,2,quantile, probs=c(0.025,0.975) )
  
  ################结果2
  HCI1=apply(Bb1,2,quantile, probs=c(0.025,0.975) )
  HCI2=apply(Bb2,2,quantile, probs=c(0.025,0.975) )
  
  ################结果3
  LHCI1=apply(Cb1,2,quantile, probs=c(0.025,0.975))
  LHCI2=apply(Cb2,2,quantile, probs=c(0.025,0.975))
  ###################
  ###################
  lower_D1=DCI1[1,];upper_D1=DCI1[2,]
  lower_D2=DCI2[1,];upper_D2=DCI2[2,]
  
  lower_L1=LCI1[1,];upper_L1=LCI1[2,]
  lower_L2=LCI2[1,];upper_L2=LCI2[2,]
  
  lower_H1=HCI1[1,];upper_H1=HCI1[2,]
  lower_H2=HCI2[1,];upper_H2=HCI2[2,]
  
  lower_LH1=LHCI1[1,];upper_LH1=LHCI1[2,]
  lower_LH2=LHCI2[1,];upper_LH2=LHCI2[2,]
  ###############################
  #########################Coverge:
  if(lower_D1[1]<B1[1]&upper_D1[1]>B1[1] )
  {Dcov11=Dcov11+1}
  if(lower_D1[2]<B1[2]&upper_D1[2]>B1[2] )
  {Dcov12=Dcov12+1}
  if(lower_D1[3]<B1[3]&upper_D1[3]>B1[3] )
  {Dcov13=Dcov13+1}
  if(lower_D1[4]<B1[4]&upper_D1[4]>B1[4] )
  {Dcov14=Dcov14+1}
  if(lower_D1[5]<B1[5]&upper_D1[5]>B1[5] )
  {Dcov15=Dcov15+1}
  if(lower_D1[6]<B1[6]&upper_D1[6]>B1[6] )
  {Dcov16=Dcov16+1}
  if(lower_D1[7]<B1[7]&upper_D1[7]>B1[7] )
  {Dcov17=Dcov17+1}
  if(lower_D1[8]<B1[8]&upper_D1[8]>B1[8] )
  {Dcov18=Dcov18+1}
  if(lower_D1[9]<B1[9]&upper_D1[9]>B1[9] )
  {Dcov19=Dcov19+1}
  if(lower_D1[10]<B1[10]&upper_D1[10]>B1[10] )
  {Dcov110=Dcov110+1}
  if(lower_D1[11]<B1[11]&upper_D1[11]>B1[11] )
  {Dcov111=Dcov111+1}
  if(lower_D1[12]<B1[12]&upper_D1[12]>B1[12] )
  {Dcov112=Dcov112+1}
  if(lower_D1[13]<B1[13]&upper_D1[13]>B1[13] )
  {Dcov113=Dcov113+1}
  if(lower_D1[14]<B1[14]&upper_D1[14]>B1[14] )
  {Dcov114=Dcov114+1}
  if(lower_D1[15]<B1[15]&upper_D1[15]>B1[15] )
  {Dcov115=Dcov115+1}
  
  #####
  if(lower_D2[1]<B2[1]&upper_D2[1]>B2[1] )
  {Dcov21=Dcov21+1}
  if(lower_D2[2]<B2[2]&upper_D2[2]>B2[2] )
  {Dcov22=Dcov22+1}
  if(lower_D2[3]<B2[3]&upper_D2[3]>B2[3] )
  {Dcov23=Dcov23+1}
  if(lower_D2[4]<rho&upper_D2[4]>rho )
  {Dcov24=Dcov24+1}
  if(lower_D2[5]<B2[5]&upper_D2[5]>B2[5] )
  {Dcov25=Dcov25+1}
  if(lower_D2[6]<B2[6]&upper_D2[6]>B2[6] )
  {Dcov26=Dcov26+1}
  if(lower_D2[7]<B2[7]&upper_D2[7]>B2[7] )
  {Dcov27=Dcov27+1}
  if(lower_D2[8]<B2[8]&upper_D2[8]>B2[8] )
  {Dcov28=Dcov28+1}
  if(lower_D2[9]<B2[9]&upper_D2[9]>B2[9] )
  {Dcov29=Dcov29+1}
  if(lower_D2[10]<B2[10]&upper_D2[10]>B2[10] )
  {Dcov210=Dcov210+1}
  if(lower_D2[11]<B2[11]&upper_D2[11]>B2[11] )
  {Dcov211=Dcov211+1}
  if(lower_D2[12]<B2[12]&upper_D2[12]>B2[12] )
  {Dcov212=Dcov212+1}
  if(lower_D2[13]<B2[13]&upper_D2[13]>B2[13] )
  {Dcov213=Dcov213+1}
  if(lower_D2[14]<B2[14]&upper_D2[14]>B2[14] )
  {Dcov214=Dcov214+1}
  if(lower_D2[15]<B2[15]&upper_D2[15]>B2[15] )
  {Dcov215=Dcov215+1}
  if(lower_D2[16]<rho&upper_D2[16]>rho )
  {Dcov216=Dcov216+1}
  
  
  ##################################
  #########################Coverge:
  if(lower_L1[1]<B1[1]&upper_L1[1]>B1[1] )
  {Lcov11=Lcov11+1}
  if(lower_L1[2]<B1[2]&upper_L1[2]>B1[2] )
  {Lcov12=Lcov12+1}
  if(lower_L1[3]<B1[3]&upper_L1[3]>B1[3] )
  {Lcov13=Lcov13+1}
  
  if(lower_L1[4]<B1[4]&upper_L1[4]>B1[4] )
  {Lcov14=Lcov14+1}
  if(lower_L1[5]<B1[5]&upper_L1[5]>B1[5] )
  {Lcov15=Lcov15+1}
  if(lower_L1[6]<B1[6]&upper_L1[6]>B1[6] )
  {Lcov16=Lcov16+1}
  
  if(lower_L1[7]<B1[7]&upper_L1[7]>B1[7] )
  {Lcov17=Lcov17+1}
  if(lower_L1[8]<B1[8]&upper_L1[8]>B1[8] )
  {Lcov18=Lcov18+1}
  if(lower_L1[9]<B1[9]&upper_L1[9]>B1[9] )
  {Lcov19=Lcov19+1}
  if(lower_L1[10]<B1[10]&upper_L1[10]>B1[10] )
  {Lcov110=Lcov110+1}
  
  if(lower_L1[11]<B1[11]&upper_L1[11]>B1[11] )
  {Lcov111=Lcov111+1}
  if(lower_L1[12]<B1[12]&upper_L1[12]>B1[12] )
  {Lcov112=Lcov112+1}
  if(lower_L1[13]<B1[13]&upper_L1[13]>B1[13] )
  {Lcov113=Lcov113+1}
  if(lower_L1[14]<B1[14]&upper_L1[14]>B1[14] )
  {Lcov114=Lcov114+1}
  if(lower_L1[15]<B1[15]&upper_L1[15]>B1[15] )
  {Lcov115=Lcov115+1}
  ###
  if(lower_L2[1]<B2[1]&upper_L2[1]>B2[1] )
  {Lcov21=Lcov21+1}
  if(lower_L2[2]<B2[2]&upper_L2[2]>B2[2] )
  {Lcov22=Lcov22+1}
  if(lower_L2[3]<B2[3]&upper_L2[3]>B2[3] )
  {Lcov23=Lcov23+1}
  
  if(lower_L2[4]<B2[4]&upper_L2[4]>B2[4] )
  {Lcov24=Lcov24+1}
  if(lower_L2[5]<B2[5]&upper_L2[5]>B2[5] )
  {Lcov25=Lcov25+1}
  if(lower_L2[6]<B2[6]&upper_L2[6]>B2[6] )
  {Lcov26=Lcov26+1}
  if(lower_L2[7]<B2[7]&upper_L2[7]>B2[7] )
  {Lcov27=Lcov27+1}
  
  if(lower_L2[8]<B2[8]&upper_L2[8]>B2[8] )
  {Lcov28=Lcov28+1}
  if(lower_L2[9]<B2[9]&upper_L2[9]>B2[9] )
  {Lcov29=Lcov29+1}
  if(lower_L2[10]<B2[10]&upper_L2[10]>B2[10] )
  {Lcov210=Lcov210+1}
  if(lower_L2[11]<B2[11]&upper_L2[11]>B2[11] )
  {Lcov211=Lcov211+1}
  
  if(lower_L2[12]<B2[12]&upper_L2[12]>B2[12] )
  {Lcov212=Lcov212+1}
  if(lower_L2[13]<B2[13]&upper_L2[13]>B2[13] )
  {Lcov213=Lcov213+1}
  if(lower_L2[14]<B2[14]&upper_L2[14]>B2[14] )
  {Lcov214=Lcov214+1}
  if(lower_L2[15]<B2[15]&upper_L2[15]>B2[15] )
  {Lcov215=Lcov215+1}
  if(lower_L2[16]<rho&upper_L2[16]>rho )
  {Lcov216=Lcov216+1}
  ####################################
  #########################Coverge:
  if(lower_H1[1]<B1[1]&upper_H1[1]>B1[1] )
  {Hcov11=Hcov11+1}
  if(lower_H1[2]<B1[2]&upper_H1[2]>B1[2] )
  {Hcov12=Hcov12+1}
  if(lower_H1[3]<B1[3]&upper_H1[3]>B1[3] )
  {Hcov13=Hcov13+1}
  
  if(lower_H1[4]<B1[4]&upper_H1[4]>B1[4] )
  {Hcov14=Hcov14+1}
  if(lower_H1[5]<B1[5]&upper_H1[5]>B1[5] )
  {Hcov15=Hcov15+1}
  if(lower_H1[6]<B1[6]&upper_H1[6]>B1[6] )
  {Hcov16=Hcov16+1}
  if(lower_H1[7]<B1[7]&upper_H1[7]>B1[7] )
  {Hcov17=Hcov17+1}
  if(lower_H1[8]<B1[8]&upper_H1[8]>B1[8] )
  {Hcov18=Hcov18+1}
  if(lower_H1[9]<B1[9]&upper_H1[9]>B1[9] )
  {Hcov19=Hcov19+1}
  if(lower_H1[10]<B1[10]&upper_H1[10]>B1[10] )
  {Hcov110=Hcov110+1}
  
  if(lower_H1[11]<B1[11]&upper_H1[11]>B1[11] )
  {Hcov111=Hcov111+1}
  if(lower_H1[12]<B1[12]&upper_H1[12]>B1[12] )
  {Hcov112=Hcov112+1}
  if(lower_H1[13]<B1[13]&upper_H1[13]>B1[13] )
  {Hcov113=Hcov113+1}
  if(lower_H1[14]<B1[14]&upper_H1[14]>B1[14] )
  {Hcov114=Hcov114+1}
  if(lower_H1[15]<B1[15]&upper_H1[15]>B1[15] )
  {Hcov115=Hcov115+1}
  ############
  if(lower_H2[1]<B2[1]&upper_H2[1]>B2[1] )
  {Hcov21=Hcov21+1}
  if(lower_H2[2]<B2[2]&upper_H2[2]>B2[2] )
  {Hcov22=Hcov22+1}
  if(lower_H2[3]<B2[3]&upper_H2[3]>B2[3] )
  {Hcov23=Hcov23+1}
  
  
  if(lower_H2[4]<B2[4]&upper_H2[4]>B2[4] )
  {Hcov24=Hcov24+1}
  if(lower_H2[5]<B2[5]&upper_H2[5]>B2[5] )
  {Hcov25=Hcov25+1}
  if(lower_H2[6]<B2[6]&upper_H2[6]>B2[6] )
  {Hcov26=Hcov26+1}
  
  if(lower_H2[7]<B2[7]&upper_H2[7]>B2[7] )
  {Hcov27=Hcov27+1}
  if(lower_H2[8]<B2[8]&upper_H2[8]>B2[8] )
  {Hcov28=Hcov28+1}
  if(lower_H2[9]<B2[9]&upper_H2[9]>B2[9] )
  {Hcov29=Hcov29+1}
  if(lower_H2[10]<B2[10]&upper_H2[10]>B2[10] )
  {Hcov210=Hcov210+1}
  if(lower_H2[11]<B2[11]&upper_H2[11]>B2[11] )
  {Hcov211=Hcov211+1}
  if(lower_H2[12]<B2[12]&upper_H2[12]>B2[12] )
  {Hcov212=Hcov212+1}
  if(lower_H2[13]<B2[13]&upper_H2[13]>B2[13] )
  {Hcov213=Hcov213+1}
  if(lower_H2[14]<B2[14]&upper_H2[14]>B2[14] )
  {Hcov214=Hcov214+1}
  if(lower_H2[15]<B2[15]&upper_H2[15]>B2[15] )
  {Hcov215=Hcov215+1}
  if(lower_H2[16]<rho&upper_H2[16]>rho )
  {Hcov216=Hcov216+1}
  
  ########################3
  if(lower_LH1[1]<B1[1]&upper_LH1[1]>B1[1] )
  {LHcov11=LHcov11+1}
  if(lower_LH1[2]<B1[2]&upper_LH1[2]>B1[2] )
  {LHcov12=LHcov12+1}
  if(lower_LH1[3]<B1[3]&upper_LH1[3]>B1[3] )
  {LHcov13=LHcov13+1}
  
  if(lower_LH1[4]<B1[4]&upper_LH1[4]>B1[4] )
  {LHcov14=LHcov14+1}
  if(lower_LH1[5]<B1[5]&upper_LH1[5]>B1[5] )
  {LHcov15=LHcov15+1}
  if(lower_LH1[6]<B1[6]&upper_LH1[6]>B1[6] )
  {LHcov16=LHcov16+1}
  
  if(lower_LH1[7]<B1[7]&upper_LH1[7]>B1[7] )
  {LHcov17=LHcov17+1}
  if(lower_LH1[8]<B1[8]&upper_LH1[8]>B1[8] )
  {LHcov18=LHcov18+1}
  if(lower_LH1[9]<B1[9]&upper_LH1[9]>B1[9] )
  {LHcov19=LHcov19+1}
  if(lower_LH1[10]<B1[10]&upper_LH1[10]>B1[10] )
  {LHcov110=LHcov110+1}
  
  if(lower_LH1[11]<B1[11]&upper_LH1[11]>B1[11] )
  {LHcov111=LHcov111+1}
  if(lower_LH1[12]<B1[12]&upper_LH1[12]>B1[12] )
  {LHcov112=LHcov112+1}
  if(lower_LH1[13]<B1[13]&upper_LH1[13]>B1[13] )
  {LHcov113=LHcov113+1}
  if(lower_LH1[14]<B1[14]&upper_LH1[14]>B1[14] )
  {LHcov114=LHcov114+1}
  if(lower_LH1[15]<B1[15]&upper_LH1[15]>B1[15] )
  {LHcov115=LHcov115+1}
  #######
  if(lower_LH2[1]<B2[1]&upper_LH2[1]>B2[1] )
  {LHcov21=LHcov21+1}
  if(lower_LH2[2]<B2[2]&upper_LH2[2]>B2[2] )
  {LHcov22=LHcov22+1}
  if(lower_LH2[3]<B2[3]&upper_LH2[3]>B2[3] )
  {LHcov23=LHcov23+1}
  
  if(lower_LH2[4]<B2[4]&upper_LH2[4]>B2[4] )
  {LHcov24=LHcov24+1}
  if(lower_LH2[5]<B2[5]&upper_LH2[5]>B2[5] )
  {LHcov25=LHcov25+1}
  if(lower_LH2[6]<B2[6]&upper_LH2[6]>B2[6] )
  {LHcov26=LHcov26+1}
  
  if(lower_LH2[7]<B2[7]&upper_LH2[7]>B2[7] )
  {LHcov27=LHcov27+1}
  if(lower_LH2[8]<B2[8]&upper_LH2[8]>B2[8] )
  {LHcov28=LHcov28+1}
  if(lower_LH2[9]<B2[9]&upper_LH2[9]>B2[9] )
  {LHcov29=LHcov29+1}
  
  if(lower_LH2[10]<B2[10]&upper_LH2[10]>B2[10] )
  {LHcov210=LHcov210+1}
  if(lower_LH2[11]<B2[11]&upper_LH2[11]>B2[11] )
  {LHcov211=LHcov211+1}
  
  if(lower_LH2[12]<B2[12]&upper_LH2[12]>B2[12] )
  {LHcov212=LHcov212+1}
  if(lower_LH2[13]<B2[13]&upper_LH2[13]>B2[13] )
  {LHcov213=LHcov213+1}
  if(lower_LH2[14]<B2[14]&upper_LH2[14]>B2[14] )
  {LHcov214=LHcov214+1}
  if(lower_LH2[15]<B2[15]&upper_LH2[15]>B2[15] )
  {LHcov215=LHcov215+1}
  if(lower_LH2[16]<rho&upper_LH2[16]>rho )
  {LHcov216=LHcov216+1}
}

#####Coverage probability of the DMC method
c(Dcov11,Dcov12,Dcov13,Dcov14,Dcov15,Dcov16, Dcov17,Dcov18,Dcov19,Dcov110,Dcov111,Dcov112,Dcov113,Dcov114,Dcov115)/NN
c(Dcov21,Dcov22,Dcov23,Dcov24,Dcov25,Dcov26,Dcov27,Dcov28,Dcov29,Dcov210,Dcov211,Dcov212,Dcov213,Dcov214,Dcov215,Dcov216)/NN
###########Coverage probability of the LSUR method
c(Lcov11,Lcov12,Lcov13,Lcov14,Lcov15,Lcov16,Lcov17,Lcov18,Lcov19,Lcov110,Lcov111,Lcov112,Lcov113,Lcov114,Lcov115)/NN
c(Lcov21,Lcov22,Lcov23,Lcov24,Lcov25, Lcov26,Lcov27,Lcov28,Lcov29,Lcov210,Lcov211, Lcov212,Lcov213,Lcov214,Lcov215,Lcov216)/NN
###########Coverage probability of the HSUR method
c(Hcov11,Hcov12,Hcov13,Hcov14,Hcov15,Hcov16,Hcov17,Hcov18,Hcov19,Hcov110, Hcov111,Hcov112,Hcov113,Hcov114,Hcov115)/NN
c(Hcov21,Hcov22,Hcov23,Hcov24,Hcov25,Hcov26,Hcov27,Hcov28,Hcov29,Hcov210,Hcov211, Hcov212,Hcov213,Hcov214,Hcov215,Hcov216)/NN
###########Coverage probability of the LHSUR method
c(LHcov11,LHcov12,LHcov13,LHcov14,LHcov15,LHcov16,LHcov17,LHcov18,LHcov19,LHcov110,LHcov111,LHcov112,LHcov113,LHcov114,LHcov115)/NN
c(LHcov21,LHcov22,LHcov23,LHcov24,LHcov25,LHcov26,LHcov27,LHcov28,LHcov29,LHcov210,LHcov211,LHcov212,LHcov213,LHcov214,LHcov215,LHcov216)/NN

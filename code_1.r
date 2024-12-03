##############
#The main  R code of paper "Bayesian Seemingly Unrelated Regression Models"
#Authors:Yang Yang, Lichun Wang and Lichun Wang
#Methods:DMC, LSUR, HSUR, LHSUR
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

n<-200#sample size
p<-3#10
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
#B1 <- c(c(0.3,0.5,0.7),rep(0, p-3)) 
#B2 <- c(c(0.2,0.4,0.6),rep(0, p-3))
B1 <- c(c(3,5,7),rep(0, p-3)) 
B2 <- c(c(2,4,6),rep(0, p-3))
#B1 <- c(c(3,3,3),rep(0, p-3)) 
#B2 <- c(c(-3,-3,-3),rep(0, p-3))
SS <-matrix(c(1,0.8,0.8,1),2,2)
er <- prnorm(rep(0,len=2),SS,n,2)
Y_1 <- X_1%*%B1+er[,1]
Y_2 <- X_2%*%B2+er[,2]
rho=SS[1,2]/sqrt(SS[1,1]*SS[2,2])
Niter=10000
## DMC method：
X <- matrix(0,ncol=length(B1)+length(B2),nrow=2*n)
X[1:n,1:length(B1)] <- X_1
X[(n+1):(2*n),(length(B1)+1):(length(B1)+length(B2))]<-X_2
y <- c(Y_1,Y_2)

hat_B1 <- solve(t(X_1)%*%X_1)%*%t(X_1)%*%Y_1
n1 <- n-p
n2 <- n-p-1
m=6000
omiga=matrix(0,m,4)
DMC_b1=matrix(0,m,length(B1))
DMC_b2=matrix(0,m,length(B2)+1)
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
  #beta_rho[i,]=c(Beta1,Beta2)
  DMC_b1[i,]=Beta1
  DMC_b2[i,]=Beta2
  #print(B)
  #print(S)
}
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)

apply(DMC_b1,2,mean)
apply(DMC_b2,2,mean)

#######################################
## Lasso SUR method：
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
  #print(Beta1.out[i,])
  Beta2.out[i,]<-updateBeta2(Y_2,X_2,Sig2.out[i-1],
                             Lambda2.out[i-1,],Y_1, X_1,Beta1.out[i,] )
  print(Beta2.out[i,length(B2)+2-1])
  Sig1.out[i]<- updateSig_1_2(Y_1,X_1,
                              Lambda1.out[i-1,], Beta1.out[i,])
  #print(Sig1.out[i])
  Sig2.out[i]<- updateSig_2_2(Y_2,X_2, Lambda2.out[i-1,],
                              Beta2.out[i,], Y_1,X_1,Beta1.out[i,])
  #print(Sig2.out[i])
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
###############################
############# HSUR  method：
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

############### LHSUR  method：
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
  return(1/rinvgauss(1,a.post,b.post))
}

updatelambda2<-function(Beta2,lambda02,Sig_2_2)
{ 
  a.post=sqrt(lambda02*Sig_2_2)/sqrt(Beta2*Beta2)
  b.post=lambda02/2
  return(1/rinvgauss(1,a.post,b.post))###
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
#  Tau.out <- matrix(0,  Niter)
Lambda01.out <- matrix(0, Niter )
Lambda02.out <- matrix(0, Niter)
################################################
Lambda1.out <- matrix(0, Niter, length(B1))
Lambda2.out <- matrix(0, Niter, length(B2))
Tau2.out <- matrix(0,  Niter)#
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

##   Gibbs sampler for LHSUR
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
  ########
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



################结果1
DCIb1=apply(DMC_b1,2,quantile, probs=c(0.025,0.975) )
DCIb2=apply(DMC_b2,2,quantile, probs=c(0.025,0.975) )

################结果2
LCIb1=apply(b1,2,quantile, probs=c(0.025,0.975) )
LCIb2=apply(b2,2,quantile, probs=c(0.025,0.975) )

################结果3
HCIb1=apply(Bb1,2,quantile, probs=c(0.025,0.975) )
HCIb2=apply(Bb2,2,quantile, probs=c(0.025,0.975) )

################结果4
LHCIb1=apply(Cb1,2,quantile, probs=c(0.025,0.975))
LHCIb2=apply(Cb2,2,quantile, probs=c(0.025,0.975))
###################

DCIb1;
DCIb2
LCIb1;
LCIb2
HCIb1;
HCIb2
LHCIb1;
LHCIb2

###############MSE
NS=Niter-Nburn
DMC_error1=DMC_b1-kronecker(matrix(B1,1),matrix(rep(1,NS),NS))
DMC_error2=DMC_b2-kronecker(matrix(c(B2,rho),1),matrix(rep(1,NS),NS))


L_error1=b1-kronecker(matrix(B1,1),matrix(rep(1,NS),NS))
L_error2=b2-kronecker(matrix(c(B2,rho),1),matrix(rep(1,NS),NS))

H_error1=Bb1-kronecker(matrix(B1,1),matrix(rep(1,NS),NS))
H_error2=Bb2-kronecker(matrix(c(B2,rho),1),matrix(rep(1,NS),NS))

LH_error1=Cb1-kronecker(matrix(B1,1),matrix(rep(1,NS),NS))
LH_error2=Cb2-kronecker(matrix(c(B2,rho),1),matrix(rep(1,NS),NS))


DMC_MSE1=diag(t(DMC_error1)%*%DMC_error1)/NS
DMC_MSE2=diag(t(DMC_error2)%*%DMC_error2)/NS
sd_DMC1=apply(DMC_b1,2,sd)
sd_DMC2=apply(DMC_b2,2,sd)

L_MSE1=diag(t(L_error1)%*%L_error1)/NS
L_MSE2=diag(t(L_error2)%*%L_error2)/NS
sd_L1=apply(b1,2,sd)
sd_L2=apply(b2,2,sd)

H_MSE1=diag(t(H_error1)%*%H_error1)/NS
H_MSE2=diag(t(H_error2)%*%H_error2)/NS
sd_H1=apply(Bb1,2,sd)
sd_H2=apply(Bb2,2,sd)

LH_MSE1=diag(t(LH_error1)%*%LH_error1)/NS
LH_MSE2=diag(t(LH_error2)%*%LH_error2)/NS
sd_LH1=apply(Cb1,2,sd)
sd_LH2=apply(Cb2,2,sd)


DMC_MSE1;L_MSE1;H_MSE1;LH_MSE1;
DMC_MSE2;L_MSE2;H_MSE2;LH_MSE2

sum(DMC_MSE1);sum(L_MSE1);sum(H_MSE1);sum(LH_MSE1);
sum(DMC_MSE2);sum(L_MSE2);sum(H_MSE2);sum(LH_MSE2);

sum(DMC_MSE1)+sum(DMC_MSE2);sum(L_MSE1)+sum(L_MSE2);sum(H_MSE1)+sum(H_MSE2);sum(LH_MSE1)+sum(LH_MSE2);
###############################3
#表1 MSE and sd
#beta1
c(sum(DMC_MSE1),sum(sd_DMC1))
c(sum(L_MSE1),sum(sd_L1))
c(sum(H_MSE1),sum(sd_H1))
c(sum(LH_MSE1),sum(sd_LH1))
#beta2
c(sum(DMC_MSE2[1:p]),sum(sd_DMC2[1:p]))
c(sum(L_MSE2[1:p]),sum(sd_L2[1:p]))
c(sum(H_MSE2[1:p]),sum(sd_H2[1:p]))
c(sum(LH_MSE2[1:p]),sum(sd_LH2[1:p]))
#rho
c(DMC_MSE2[1+p],sd_DMC2[1+p])
c(L_MSE2[1+p],sd_L2[1+p])
c(H_MSE2[1+p],sd_H2[1+p])
c(LH_MSE2[1+p],sd_LH2[1+p])
#################################
#beta1  mse sd
print(matrix(c(sum(DMC_MSE1),sum(sd_DMC1),
               sum(L_MSE1),sum(sd_L1),
               sum(H_MSE1),sum(sd_H1),
               sum(LH_MSE1),sum(sd_LH1)),4,2,2))
#beta2  mse sd
print(matrix(c(sum(DMC_MSE2[1:p]),sum(sd_DMC2[1:p]),
               sum(L_MSE2[1:p]),sum(sd_L2[1:p]),
               sum(H_MSE2[1:p]),sum(sd_H2[1:p]),
               sum(LH_MSE2[1:p]),sum(sd_LH2[1:p])),4,2,2))
#rho mse sd
print(matrix(c(DMC_MSE2[1+p],sd_DMC2[1+p],
               L_MSE2[1+p],sd_L2[1+p],
               H_MSE2[1+p],sd_H2[1+p],
               LH_MSE2[1+p],sd_LH2[1+p]),4,2,2))
#sum
sum(DMC_MSE1)+sum(DMC_MSE2);sum(L_MSE1)+sum(L_MSE2);sum(H_MSE1)+sum(H_MSE2);sum(LH_MSE1)+sum(LH_MSE2);

# Geweke's test for LSUR
geweke.diag(mcmc(b1))$z
geweke.diag(mcmc(b2))$z
# Geweke's test for HSUR
geweke.diag(mcmc(Bb1))$z
geweke.diag(mcmc(Bb2))$z
# Geweke's test for LHSUR
geweke.diag(mcmc(Cb1))$z
geweke.diag(mcmc(Cb2))$z


#beta1  mse sd
print(matrix(c(sum(DMC_MSE1),sum(sd_DMC1),
               sum(L_MSE1),sum(sd_L1),
               sum(H_MSE1),sum(sd_H1),
               sum(LH_MSE1),sum(sd_LH1),
               sum(DMC_MSE2[1:p]),sum(sd_DMC2[1:p]),
               sum(L_MSE2[1:p]),sum(sd_L2[1:p]),
               sum(H_MSE2[1:p]),sum(sd_H2[1:p]),
               sum(LH_MSE2[1:p]),sum(sd_LH2[1:p]),
               DMC_MSE2[1+p],sd_DMC2[1+p],
               L_MSE2[1+p],sd_L2[1+p],
               H_MSE2[1+p],sd_H2[1+p],
               LH_MSE2[1+p],sd_LH2[1+p]),12,2,2))
########### MSE(b)
sum(DMC_MSE1)+sum(DMC_MSE2);sum(L_MSE1)+sum(L_MSE2);sum(H_MSE1)+sum(H_MSE2);sum(LH_MSE1)+sum(LH_MSE2);






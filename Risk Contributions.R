#Initialization for the proposal value.
library(MASS)
library(LaplacesDemon)
library(mcmcse)
library(dplyr)

#This will house the results
estimators<-as.data.frame(matrix(nrow=3, ncol=4))
colnames(estimators)<-c("MC", "NW", "GR", "MH")

#Covariance matrix
sig<-matrix(c(1, -0.5, 0.3, -0.5, 1, 0.5, 0.3, 0.5, 1), nrow=3)


####Generate 10,000 values from the multivariate t distribution####
x.mc<-rmvt(n=10^6, df=4, S=sig, mu=rep(0,3))
x.mc<-cbind(x.mc, apply(x.mc, 1, sum))

#Estimate VaR of entire portfolio, at p=0.999, using the method on p.1585 - although we want right tail
v<-x.mc[order(x.mc[,4], decreasing = T)[1000], 4]

#### MC Estimator####

start<-Sys.time()

#We want about 100 observations in the acceptance interval
del.int<-x.mc[order(x.mc[,4], decreasing = T)[500:1500], 4]

#Find the observations such that S is in the acceptance interval
mc.sub<-x.mc[x.mc[,4] %in% del.int, 1:4]

#Save the estimators
estimators[,1]<-apply(mc.sub[, 1:3], 2, mean)

Sys.time()-start

#mat<-(t(mc.sub[,1:3]-mc.mat) %*% (mc.sub[,1:3]-mc.mat))

#### NW Estimator####
start<-Sys.time()
#Estimate the kernel (see p. 1586-1587)
alt_j<-1.06*sd(x.mc[,4])*10^(-0.2*6)
kers<-dnorm((x.mc[,4]-v)/alt_j)

#Save the estimators
estimators[,2]<-apply(x.mc[,1:3], 2, function(x) sum(x*kers)/sum(kers))

Sys.time()-start


#### GR Method ####

start<-Sys.time()

for(i in 1:3){
  mean.x<-mean(x.mc[,i])
  mean.s<-mean(x.mc[,4])
  
  #Calculate coefficients
  b1<-sum((x.mc[,i]-mean.x)*(x.mc[,4]-mean.s))/sum((x.mc[,4]-mean.s)^2)
  b0<-mean.x-b1*mean.s
  
  #Return estimators
  estimators[i,3]<-b0+b1*v
}
Sys.time()-start

#### MH Method ####

start<-Sys.time()

#Parameters
mc.mat<-t(as.matrix(estimators[,1]))
mc.mat<-mc.mat[rep(1, 1001),]
sig.mh<-(t(mc.sub[,1:3]-mc.mat) %*% (mc.sub[,1:3]-mc.mat))/nrow(mc.sub)
mu<-c(estimators[1:2,3], v-estimators[1,3]-estimators[2,3])

#This function will be used for the MpCN distribution
norm.mpcn<-function(x, sig){
  sig.eigen<-eigen(sig)
  sig.sqrt<-sig.eigen$vectors%*%diag(sqrt(sig.eigen$values))%*%t(sig.eigen$vectors)
  (norm(solve(sig.sqrt)%*%(x-mu))^2)/2
}

#Generate a value from the MpCN distribution
mpcn<-function(mu, p, x, sig){
  
  #Generate value from W
  w<-rmvn(n=1, mu = rep(0, 3), S = sig)
  
  #Generate value from Z
  s.par<-norm.mpcn(x,sig)
  z<-rgamma(n=1, shape = length(x)/2, scale = s.par)
  
  mu+sqrt(p)*(x-mu)+sqrt(1-p)*sqrt(z)^(-1)*w
}

#This matrix will house the observations
x.mh<-matrix(nrow=10^6, ncol=3)

#First value
x.can<-mpcn( mu, 0.8, rep(as.numeric(v/3), 3), sig.mh)

#Components of the MH algorithm (see p.1584)
p.new<-dmvt(x.can, mu, S=sig.mh, df=4)
p.old<-dmvt(rep(as.numeric(v/3), 3), mu, S=sig.mh, df=4)
q.new<-norm.mpcn(as.numeric(x.can), sig.mh)
q.old<-norm.mpcn(rep(as.numeric(v/3), 3), sig.mh)

#MH Algorithm
candidate<-min(1, (p.new/p.old)*(q.old/q.new)^(-3))
if(candidate>runif(1)){
  x.mh[1,]<-x.can
}else{x.mh[1,]<-rep(v/3, 3)}

#Other 9999 values
for(i in 2:(10^6)){
  x.can<-mpcn(x = x.mh[i-1,], mu = mu, sig = sig.mh, p=0.8)
  
  p.new<-dmvt(x.can, mu, S=sig.mh, df=4)
  p.old<-dmvt(x.mh[i-1,], mu, S=sig.mh, df=4)
  q.new<-norm.mpcn(as.numeric(x.can), sig.mh)
  q.old<-norm.mpcn(x.mh[i-1,], sig.mh)
  
  candidate<-min(1, (p.new/p.old)*(q.old/q.new)^(-3))
  if(candidate>runif(1)){
    x.mh[i,]<-x.can
  }else{x.mh[i,]<-x.mh[i-1,]}
}

estimators[,4]<-colMeans(x.mh)

Sys.time()-start


####Standard Errors####

##Matrix to house the observations
SE<-as.data.frame(matrix(ncol=3, nrow=3))
colnames(SE)<-c("MC", "GR", "MH")

#Crude Monte Carlo
mc.mat<-t(as.matrix(estimators[,1]))
mc.mat<-mc.mat[rep(1, 1001),]
sig.mh<-(t(mc.sub[,1:3]-mc.mat) %*% (mc.sub[,1:3]-mc.mat))/nrow(mc.sub)
SE[,1]<-sqrt(diag(sig.mh)/nrow(mc.sub))

#Generalized Regression
GR.Res<-matrix(ncol=3, nrow=nrow(x.mc))

Y<-cbind(rep(1, nrow(x.mc)), x.mc[,4])
Sig.gr.constant<-((t(Y) %*% Y)/nrow(x.mc)) %>% solve

for(i in 1:3){
  mean.x<-mean(x.mc[,i])
  mean.s<-mean(x.mc[,4])
  b1<-sum((x.mc[,i]-mean.x)*(x.mc[,4]-mean.s))/sum((x.mc[,4]-mean.s)^2)
  b0<-mean.x-b1*mean.s
  
  GR.Res<-x.mc[,i]-b0-b1*-x.mc[,4]
  
  Sig.gr<-sd(GR.Res)*Sig.gr.constant
  
  SE[i,2]<-sqrt(Sig.gr[1,1]+2*v*Sig.gr[1,2]+v^2*Sig.gr[2,2])/sqrt(nrow(x.mc))
}

#MCMC
SE[,3]<-mcse.mat(x.mh, method = "bm")[,2]

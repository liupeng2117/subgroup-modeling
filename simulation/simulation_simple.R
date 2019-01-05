#Simulation
#============================
# 500 subjects, 2 subgroups, 2 clinical variables x1, x2, 2 genes G1, G2
#============================
set.seed(123)
n.rep=100
par.all<-matrix(,ncol=11,nrow=n.rep)
for(rep in 1:n.rep){
n=500 # sample size
T=c(rep(0,n/2),rep(1,n/2)) #treatment
X1=rnorm(n,2,1) # clinical var1
X2=rnorm(n,4,2) # clinical var2
G1=rnorm(n,0,1) # gene 1
G2=rnorm(n,1,2)

alpha=c(1, 3) # coefficient of treatment
beta1=c(0.4,0.6) # coefficient of clinical var1
beta2=c(0.8,1.0) # coefficient of clinical var2
gamma1=1.2 # coefficient of gene 1
gamma2=0.5 # coefficient of gene 2
miu=1 # overall mean
sigma2=1 # overall variance

pai=exp(G1*gamma1+G2*gamma2)/(1+exp(G1*gamma1+G2*gamma2))
Z=sapply(1:length(pai), function(i) rbinom(1,size=1,prob=1-pai[i]))
Z=Z+1

e=rnorm(n,0,sigma2)

#calculate outcome by regression function
Y=miu+alpha[Z]*T+beta1[Z]*X1+beta2[Z]*X2+e

#E-M algorithm

#10 initials and choose the best
n.int=10
est.par=matrix(,ncol=11,nrow=n.int)
for(i in 1:n.int){
#initial values
alpha_int=runif(2,0,3)
beta1_int=runif(2,0,3)
beta2_int=runif(2,0,3)
gamma1_int=runif(1,0,1)
gamma2_int=runif(1,0,1)
miu_int=runif(1,0,1)
sigma2_int=runif(1,1,3)
theta=c(alpha_int, beta1_int, beta2_int, gamma1_int, gamma2_int,miu_int, sigma2_int)

pai_int=exp(G1*gamma1_int+G2*gamma2_int)/(1+exp(G1*gamma1_int+G2*gamma2_int))
f1_int=(1/sqrt(2*pi*sigma2_int))*exp(-(Y-miu_int-T*alpha_int[1]-X1*beta1_int[1]-X2*beta2_int[1])^2/(2*sigma2_int))
f2_int=(1/sqrt(2*pi*sigma2_int))*exp(-(Y-miu_int-T*alpha_int[2]-X1*beta1_int[2]-X2*beta2_int[2])^2/(2*sigma2_int))
(ll_int=sum(log(pai_int* f1_int+(1-pai_int)*f2_int)))

l=1
repeat{
  alpha_old=theta[1:2]
  beta1_old=theta[3:4]
  beta2_old=theta[5:6]
  gamma1_old=theta[7]
  gamma2_old=theta[8]
  miu_old=theta[9]
  sigma2_old=theta[10]
  #==E-STEP==#
  pai_old=exp(G1*gamma1_old+G2*gamma2_old)/(1+exp(G1*gamma1_old+G2*gamma2_old))
  f1_old=(1/sqrt(2*pi*sigma2_old))*exp(-(Y-miu_old-T*alpha_old[1]-X1*beta1_old[1]-X2*beta2_old[1])^2/(2*sigma2_old))
  f2_old=(1/sqrt(2*pi*sigma2_old))*exp(-(Y-miu_old-T*alpha_old[2]-X1*beta1_old[2]-X2*beta2_old[2])^2/(2*sigma2_old))
  
  #calculate the expected value of Z
  w_old=(pai_old*f1_old)/(pai_old*f1_old+(1-pai_old)*f2_old)
  #==M-STEP==#
  #update gamma1 gamma2 by Newton-raphson algorithm 
  func <- function(x,G,gamma,w,i) {
    sum(G[,i]*(w-exp(G[,i]*x+G[,-i]*gamma[-i])/(1+exp(G[,i]*x+G[,-i]*gamma[-i]))))
  }
  gamma1_new<-uniroot(func,interval=c(-10,10), G=cbind(G1,G2),gamma=c(gamma1_old, gamma2_old),w=w_old,i=1)$root
  gamma2_new<-uniroot(func,interval=c(-10,10), G=cbind(G1,G2),gamma=c(gamma1_new, gamma2_old),w=w_old,i=2)$root
  #update miu, alpha
  miu_new=sum(w_old*(Y-X1*beta1_old[1]-X2*beta2_old[1]-T*alpha_old[1])+
                   (1-w_old)*(Y-X1*beta1_old[2]-X2*beta2_old[2]-T*alpha_old[2]))/n
  alpha_new=alpha_old
  alpha_new[1]<- sum(w_old*T*(Y-miu_new-X1*beta1_old[1]-X2*beta2_old[1]))/sum(w_old*T^2)
  alpha_new[2]<- sum((1-w_old)*T*(Y-miu_new-X1*beta1_old[2]-X2*beta2_old[2]))/sum((1-w_old)*T^2)
  #update beta1 beta2
  beta1_new=beta1_old
  beta1_new[1]<-sum(w_old*X1*(Y-miu_new-T*alpha_new[1]-X2*beta2_old[1]))/sum(w_old*X1^2)
  beta1_new[2]<-sum((1-w_old)*X1*(Y-miu_new-T*alpha_new[2]-X2*beta2_old[2]))/sum((1-w_old)*X1^2)
  beta2_new=beta2_old
  beta2_new[1]<-sum(w_old*X2*(Y-miu_new-T*alpha_new[1]-X1*beta1_new[1]))/sum(w_old*X2^2)
  beta2_new[2]<-sum((1-w_old)*X2*(Y-miu_new-T*alpha_new[2]-X1*beta1_new[2]))/sum((1-w_old)*X2^2)

  #update sigma2
  sigma2_new=sum(w_old*(Y-miu_new-T*alpha_new[1]-X1*beta1_new[1]-X2*beta2_new[1])^2+
                 (1-w_old)*(Y-miu_new-T*alpha_new[2]-X1*beta1_new[2]-X2*beta2_new[2])^2)/n
  
  theta_new=c(alpha_new, beta1_new, beta2_new, gamma1_new, gamma2_new,miu_new, sigma2_new)
  dis=sqrt(sum((theta_new-theta)^2))
  theta=theta_new
  l=l+1
  if(dis<1*10^(-7) |l>1000){
    break
  }
}

(alpha_est=theta[1:2])
(beta1_est=theta[3:4])
(beta2_est=theta[5:6])
(gamma1_est=theta[7])
(gamma2_est=theta[8])
(miu_est=theta[9])
(sigma2_est=theta[10])

pai_est=exp(G1*gamma1_est+G2*gamma2_est)/(1+exp(G1*gamma1_est+G2*gamma2_est))
f1_est=(1/sqrt(2*pi*sigma2_est))*exp(-(Y-miu_est-T*alpha_est[1]-X1*beta1_est-X2*beta2_est[1])/(2*sigma2_est))
f2_est=(1/sqrt(2*pi*sigma2_est))*exp(-(Y-miu_est-T*alpha_est[2]-X1*beta1_est-X2*beta2_est[2])/(2*sigma2_est))
(ll=log(prod(pai_est* f1_est+(1-pai_est)*f2_est)))

est.par[i,]<-c(theta,ll)
}
#choose the best par estimates which has the maximum ll
best.par<-est.par[which.max(est.par[,11]),]
par.all[rep,]<-best.par
print(rep)
print(best.par)
}
for(i in 1:n.rep){
  if(par.all[i,1]>par.all[i,2]) {
    par.all[i,]<-c(par.all[i,c(2,1,4,3,6,5)], -par.all[i,7:8], par.all[i,c(9:11)])
  }
}
# mean and se of the par est
apply(par.all,2,mean)
apply(par.all,2,sd)

pai.all=exp(G1 %o% par.all[,7] + G2 %o% par.all[,8])/(1+exp(G1 %o% par.all[,7]+G2 %o% par.all[,8]))

mean(sapply(1:n.rep, function(i) cor(pai,pai.all[,i])))
plot(pai,pai.all[,1],xlab="true pai_i", ylab="estimated pai_i", pch=16,cex=0.5)
abline(0,1,col=2)





#============================
func2 <- function(x,G,gamma,i) {
  sum(G[,i]*(0.4-exp(G[,i]*x+G[,-i]*gamma[-i])/(1+exp(G[,i]*x+G[,-i]*gamma[-i]))))
}
G=matrix(c(1,2,2,1,1,2),ncol=2)
gamma=c(0.3,0.5)
i=1
func3<- function(x,G=G,gamma=gamma,i=1) {
  Reduce("+", lapply(1:3,function(j) G[j,i]*(0.4-exp(G[j,i]*x+G[j,-i]*gamma[-i])/(1+exp(G[j,i]*x+G[j,-i]*gamma[-i])))))
}
curve(func3,xlim=c(-3,3), col='blue', lwd=2, lty=2, ylab='f(x)')
uniroot(func2,interval=c(-3,3), G=G,gamma=gamma,i=1)
library(rootSolve)
uniroot.all(function(x) func3(x,G=G,gamma=gamma,i=1),interval=c(-3, 10))

#library(numDeriv)
#genD(func = func2, x = 2)
pai=exp(G1*gamma1+G2*gamma2)/(1+exp(G1*gamma1+G2*gamma2))
f1=(1/sqrt(2*pi*sigma2))*exp(-(Y-miu-T*alpha[1]-X1*beta1[1]-X2*beta2[1])/(2*sigma2))
f2=(1/sqrt(2*pi*sigma2))*exp(-(Y-miu-T*alpha[2]-X1*beta1[2]-X2*beta2[2])/(2*sigma2))
(ll=log(prod(pai* f1+(1-pai)*f2)))
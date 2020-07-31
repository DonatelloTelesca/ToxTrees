setwd(".../tree_splines_biostatistics/project")

library(splines)
library(MASS)

data = read.csv("test_simul_17.txt",header=FALSE,sep=" ")
# data = read.csv("test_simul_25.txt",header=FALSE,sep=" ")
# data = read.csv("test_simul_32.txt",header=FALSE,sep=" ")
# data = read.csv("test_simul_36.txt",header=FALSE,sep=" ")

beta = data[1:1999,]
beta = apply(beta,2,as.character)
beta = apply(beta,2,as.numeric)

tau2 = data[2009:22007,1]
tau2 = as.numeric(as.character(tau2))

# 22008+100000/13
phid = data[22008:29700,]
phid = as.vector(t(phid))
phid = phid[1:100000]
phid = as.numeric(as.character(phid))

sigma2 = data[29701:49699,1]
tmp = as.character(sigma2)
tmp2 = rep("",19999)
for (i in 1:19999)
{
  tmp2[i] = unlist(strsplit(tmp[i],"(",fixed=TRUE))[[3]]
}

tmp3 = rep("",19999)
for (i in 1:19999)
{
  tmp3[i] = unlist(strsplit(tmp2[i],")",fixed=TRUE))[[1]]
}

sigma2 = as.numeric(tmp3)

B = bs(x=seq(1,11,by=1),knots = c(2:10), degree = 3, intercept = TRUE, Boundary.knots = c(1,11))
dmax=11
post = matrix(0,nrow=1999,ncol=11)

for (it in 1:1999)
{
  phid.tmp = phid[80001+(it)*10]
  sigma2.tmp = sigma2[(it)*10+1]
  sig.tmp = diag(dmax)
  sig.tmp <- sigma2.tmp * phid.tmp^abs(row(sig.tmp)-col(sig.tmp))
  beta.tmp = beta[it,]
  post[it,] = mvrnorm(1,B %*% beta.tmp,sig.tmp)
}

data.curves = read.csv("simul_dose2.csv",sep=" ", header=FALSE)
dose = c(0,0.39,0.78,1.56,3.125,6.25,12.5,25,50,100,200)

y.mean = apply(post,2,mean)
y.5 = apply(post,2,quantile,probs=0.05)
y.95 = apply(post,2,quantile,probs=0.95)

i = 17
# i = 25
# i = 32
# i = 36

plot(1,1,col="white",ylim=range(data.curves),xlim=c(0.39,200),xlab="Dose",ylab="Response", main=i,log="x")
points(dose,y.mean,type="l",lwd=3)
points(dose,y.5,type="l",lwd=3)
points(dose,y.95,type="l",lwd=3)
points(dose,data.curves[(((i-1)*dmax)+1):(((i-1)*dmax)+dmax),1],type="p",col=i,cex=1,pch=i%%26)
points(dose,data.curves[(((i-1)*dmax)+1):(((i-1)*dmax)+dmax),2],type="p",col=i,cex=1,pch=i%%26)
points(dose, apply(data.curves[(((i-1)*dmax)+1):(((i-1)*dmax)+dmax),c(1:2)],1,mean),type="l",col=i,lwd=1)

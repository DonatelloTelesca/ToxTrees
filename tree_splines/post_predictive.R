setwd(".../tree_splines_biostatistics/project")

library(splines)
library(MASS)

# read results file
data = read.csv("test_simul_dose2_posterior_predictive.txt",header=FALSE,sep=" ")

setwd(".../tree_splines_biostatistics/data")
# read original data file
data2 = read.csv("simul_dose2.csv",sep=" ", header=FALSE)

# parameters
N = 20000
M = 100000
n_step = 10
I = 50
P = 8
K = 2
dmax = 11
delta = 2

# splines coefficients
end.beta = ((N/n_step)-1)*I

beta = data[1:end.beta,]
beta = apply(beta,2,as.character)
beta = apply(beta,2,as.numeric)

# smoothing parameter for splines
begin.tau2 = end.beta+P+2
end.tau2 = begin.tau2 + N -2

tau2 = data[begin.tau2:end.tau2,1]
tau2 = as.numeric(as.character(tau2))

plot(tau2,type="l",ylab="tau2")

# dose correlation
begin.phid = end.tau2+1
end.phid = floor(begin.phid + M/(dmax+delta))

phid = data[begin.phid:end.phid,]
phid = as.vector(t(phid))
phid = phid[1:M]
phid = as.numeric(as.character(phid))

plot(phid,type="l",ylab="phid")

# variance
begin.sigma2 = end.phid+1
end.sigma2 = begin.sigma2+N-2
sigma2 = data[begin.sigma2:end.sigma2,1]
sigma2 = as.character(sigma2)

for (i in 1:N-1)
{
  sigma2[i] = unlist(strsplit(sigma2[i],"(",fixed=TRUE))[[3]]
}

for (i in 1:N-1)
{
  sigma2[i] = unlist(strsplit(sigma2[i],")",fixed=TRUE))[[1]]
}

sigma2 = as.numeric(sigma2)

plot(sigma2,type="l",ylab="sigma2")

# posterior predictive curves
B = bs(x=seq(1,dmax,by=1),knots = c(2:(dmax-1)), degree = 3, intercept = TRUE, Boundary.knots = c(1,dmax))
post = matrix(0,nrow=(N/n_step-1)*I,ncol=dmax)

# generate y
for (it in 1:((N/n_step)-1))
{
  phid.tmp = phid[M-N+1+(it)*n_step]
  sigma2.tmp = sigma2[(it)*n_step+1]
  sig.tmp = diag(dmax)
  sig.tmp <- sigma2.tmp * phid.tmp^abs(row(sig.tmp)-col(sig.tmp))
  for (i in 1:I)
  {
    beta.tmp = beta[(it-1)*I+i,]
    post[(it-1)*I+i,] = mvrnorm(1,B %*% beta.tmp,sig.tmp)
  }
}

# plot curves
v = c(1:((N/n_step-1)*I))
par(mfrow=c(3,2))
for (i in 1:I)
{
  index = which(v[v %%(I+1)]==i)
  y.mean = apply(post[index,],2,mean)
  y.5 = apply(post[index,],2,quantile,probs=0.05)
  y.95 = apply(post[index,],2,quantile,probs=0.95)
  plot(1,1,col="white",ylim=range(data2),xlim=c(1,dmax),xlab="Dose",ylab=i)
  points(y.mean,type="l",lwd=3)
  points(y.5,type="l",lwd=3)
  points(y.95,type="l",lwd=3)
  for (k in 1:K)
  {
    points(data2[(((i-1)*dmax)+1):(((i-1)*dmax)+dmax),k],type="p",col=i, pch=i%%26)
  }
  points(colMeans(t(data2[(((i-1)*dmax)+1):(((i-1)*dmax)+dmax),(1:K)])),type="l",col=i)
}


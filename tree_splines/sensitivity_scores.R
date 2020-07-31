setwd(".../tree_splines_biostatistics/project")
library("splines")

beta = read.csv("test_simul_hypercube.txt",header=FALSE,sep=" ",nrows=19990000)
beta = apply(beta,2,as.character)
beta = apply(beta,2,as.numeric)

# from beta coefficients to y
B = bs(x=seq(1,11,by=1),knots = c(2:10), degree = 3, intercept = TRUE, Boundary.knots = c(1,11))
fit.spline = function(v)
{
 fit = B %*% v
 return(fit)
}

# inner product
inner = function(v)
{
  return(v%*%v)
}

# number of iterations
nb_it = 1999
# number of rows of M, M' and Mj combined
m = 10000
# number of doses
D = 11

# expectation
M_exp = matrix(0,ncol=D,nrow=nb_it)
M_exp2 = matrix(0,ncol=D,nrow=nb_it)


for (i in 1:nb_it)
{
  tmp = beta[((i-1)*m+1):((i-1)*m+1000),]
  tmp = apply(tmp,1,fit.spline)
  M_exp[i,] = apply(t(tmp),2,mean)
  M_exp2[i,] = apply(t(tmp),2,inner)/1000
}

# variance
M_var = M_exp2 -M_exp^2

# Sj and Tj
M_SM = matrix(0,ncol=D,nrow=nb_it*8)
M_SMprime = matrix(0,ncol=D,nrow=nb_it*8)

for (j in 1:8)
{
  for (i in 1:nb_it)
  {
    tmp = beta[((i-1)*m+1):((i-1)*m+1000),]
    tmp = apply(tmp,1,fit.spline)
    tmpprime = beta[((i-1)*m+1001):((i-1)*m+2000),]
    tmpprime = apply(tmpprime,1,fit.spline)
    tmpj = beta[((i-1)*m+(j+1)*1000+1):((i-1)*m+(j+2)*1000),]
    tmpj = apply(tmpj,1,fit.spline)
    M_SM[(j-1)*nb_it+i,] = diag(tmp %*% t(tmpj))/1000
    M_SMprime[(j-1)*nb_it+i,] = diag(tmpprime %*% t(tmpj))/1000
  }
}

# repeat expectation and variance for each predictor
M_exp8 = rbind(M_exp, M_exp, M_exp, M_exp, M_exp, M_exp, M_exp, M_exp)
M_var8 = rbind(M_var, M_var, M_var, M_var, M_var, M_var, M_var, M_var)

# calculate scores for each predictor
M_S = (M_SM - M_exp8^2) / M_var8
v_S = apply(M_S,1,mean)
M_T = (M_SMprime - M_exp8^2)/ M_var8
v_T = apply(M_T,1,mean)

# organize scores for plot
S_tmp = matrix(0,nrow=8,ncol=1999)
T_tmp = matrix(0,nrow=8,ncol=1999)
for (j in 1:8)
{
  S_tmp[j,] = abs(v_T[(nb_it*(j-1)+1):(nb_it*j),1])
  T_tmp[j,] = abs(1-v_S[(nb_it*(j-1)+1):(nb_it*j),1])
}

# plot scores
par(mfrow=c(2,1))
boxplot(t(S_tmp),ylim=c(0,1),ylab="S",names=name.var,cex.axis=0.8)
boxplot(t(T_tmp),ylim=c(0,1),ylab="T",names=name.var,cex.axis=0.8)



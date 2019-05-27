install.packages("ppcor") # package needed to do partial correlation
library("ppcor")

# network
# DESCRIPTION: 
# Computes a 0/1/-1 matrix based on three options. Should centralisation be performed within sub-experiments, Should the correlation be computed with Pearson or partial correlation, should bootstrap replicated be calculated.

# USAGE: 
#  network(x, y, q=0.005, B=50, method.cent="NonCSE", method.cor="pearson", method.precision="cutoff")

# REQUIRED ARGUMENTS: 
#	x			            matrix with n samples and p variables. Missing values (NAs) are not allowed. 
#	y			            integer vector of length n, where samples having the same numbers are defined to be in the same sub-experiments.  Missing values (NAs) are allowed.
#	q			            desired fraction of edges
#	B			            number of bootstrap replicates (50 default)
#	method.cent	    	a character string indicating centralisation within sub-experiments is to be computed. One of "NonCSE" (default), or "CSE". 
#	method.cor		    a character string indicating which correlation method is to be computed. One of "pearson" (default), or "partial". 
#	method.precision	a character string indicating if the sparsity is to be computed on bootstrap replicates or not. One of "cutoff" (default), or "bootstrap". 	

#DETAILS:
# Partial correlation is the correlation of two variables while controlling for a third or more other variables. When the determinant of variance-covariance matrix is numerically zero, the Moore-Penrose generalized matrix inverse is used. 
#	NAs are not allowed 

#VALUE: 
#  Returns a matrix with specified sparsity (n rows and n columns). 

network <-function(x, y, q=0.005, B=50, method.cent="NonCSE", method.cor="pearson", method.precision="cutoff")
{
  if(method.cent=="CSE")
  {
    x <- CSE(x,y)
  }
  x <- corMatrix(x,method.cor)
  if (method.precision=="bootstrap")
  {
    x <- precision.bootstrap(x, B, q, method.cor)
  }
  if (method.precision=="cutoff")
  {
    x <- precision.cutoff(x,q)
  }
  return(x)
}

# method.cor
# DESCRIPTION: 
# Calculates Pearson or partial correlation. Partial correlation calculates the pairwise correlation for each pair of variables given others. Pearson calculates the pairwise correlation. 

# USAGE: 
# cor.matrix(x,method.cor="pearson")

# REQUIRED ARGUMENTS: 
# x			      matrix with n samples and p variables.
# method.cor	a character string indicating which correlation method is to be computed. One of "pearson" (default), or "partial". 

# DETAILS:
# Partial correlation is the correlation of two variables while controlling for a third or more other variables. When the determinant of variance-covariance matrix is numerically zero, the Moore-Penrose generalized matrix inverse is used. 
# Nas are not allowed

# VALUE: 
# Returns a matrix with specified correlation (n rows and n columns). 

corMatrix<- function(x, method.cor="pearson")
{
  if (method.cor=="partial")
  {
    x.cor <- pcor(x)$estimate 
  }
  if (method.cor=="pearson")
  {
    x.cor <- cor(x)
  }
  return(x.cor)
}

# precision.cutoff
# DESCRIPTION: 
#  Computes a 0/1/-1 matrix from a correlation matrix. The pre-defined sparsity sets the amount of non-zero edges in the matrix. For example, .05 sparsity means 5% of the edges are non-zero. The highest absolute value of the correlation matrix are set to 1/-1 depending on the correlation sign. 

# USAGE: 
# Presision.cutoff(x,q)

# REQUIRED ARGUMENTS: 
# x correlation matrix. Missing values (NAs) are not allowed. 
# q desired fraction of edges.

# DETAILS:
# NAs are not allowed. 

# VALUE: 
# Returns a matrix with specified sparcity (n rows and n columns). 

precision.cutoff <- function(x, q)
{
  n <- nrow(x)
  i=1
  while(((sum(abs(x)>(i-0.5)/1000)- n)/(n^2-n))>q)
  {
    i <- i+0.1		
  }
  x <- (abs(x)>((i-0.5)/1000))*sign(x)
  return(x)
}

# precicion.bootstrap
# DESCRIPTION: 
# Computes a 0/1/-1 matrix based on B bootstrap replicates from a correlation matrix. The pre-defined sparsity sets the amount of non-zero edges in the matrix. The highest absolute value of the correlation matrix are set to 1/-1. For example, .05 sparsity means 5% of the edges are non-zero. Calculates pearson or partial correlation. Partial correlation calculates the pairwise correlation for each pair of variables given others. Pearson calculates the the pairwise correlation. 

# USAGE: 
#  precision.bootstrap(x, B=50, q=0.005,method.cor="pearson")

# REQUIRED ARGUMENTS: 
#	x	          matrix with n samples and p variables. Missing values (NAs) are not allowed
#	B	          number of bootstrap replicates (50 default)
#	q	          desired fraction of edges
# method.cor	a character string indicating which correlation method is to be computed. One of "pearson" (default), or "partial". 

# DETAILS:
# Partial correlation is the correlation of two variables while controlling for a third or more other variables. When the determinant of variance-covariance matrix is numerically zero, the Moore-Penrose generalized matrix inverse is used. 
#	NAs are not allowed 

# VALUE: 
#  Returns a matrix with specified sparsity (n rows and n columns). 

precision.bootstrap <- function(x, B=50, q=0.005, method.cor="pearson")
{
  ag.prec <- matrix(0, nrow=ncol(x), ncol=ncol(x))      
  rownames(ag.prec)<-colnames(x) 
  colnames(ag.prec)<- colnames(x) 
  n <- length(x[,1])
  for(b in 1:B) 
  {
    BootstrapReplicates <- sample(n, replace = T)
    x.cor <- corMatrix(x[BootstrapReplicates,], method.cor)
    x.prec <- precision.cutoff(x.cor,q)
    ag.prec <- ag.prec + x.prec
  }
  x.prec <- precision.cutoff(ag.prec/B,q)
  return(x.prec)
}

# CSE
# DESCRIPTION: 
# Centralisation within sub-experiments

# USAGE: 
# CSE(x, y)

# REQUIRED ARGUMENTS: 
# x matrix with n samples and p variables. 
# y integer vector of length n, where samples having the same numbers are defined to be in the same sub-experiments.  Missing values (NAs) are allowed. 

# DETAILS:
# Sub-experiments with less than two samples will generate NAs. 

# VALUE: 
# centralisation returns a matrix with centralized data (n rows and p columns). 

# EXAMPLES: 
# x <- matrix(50:1, ncol = 10)
# y <- c(1, 3, 3, 1, 7)
# CSE(x, y)

# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] 
# [1,]  1.5  1.5  1.5  1.5  1.5  1.5  1.5  1.5  1.5   1.5
# [2,]  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5  0.5   0.5
# [3,] -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5  -0.5
# [4,] -1.5 -1.5 -1.5 -1.5 -1.5 -1.5 -1.5 -1.5 -1.5  -1.5
# [5,]   NA   NA   NA   NA   NA   NA   NA   NA   NA    NA

CSE <- function(x, y)
{
  if(sum(data.frame(table(y))[,2]==1)>0)
  {
    print("one or more subexperiments consist of one replicate. Control your y values if all sub-experiments should consist of more than 1 replicate. Otherwise, remove the the sub-experiments with only one replicate to preform centralisation")
  }
  sub.exp <- unique(y)
  for(i in 1:length(sub.exp))
  {
    if(length(x[y==sub.exp[i],1])>1)
    {
      x[y==sub.exp[i],] <- t(t(x[y==sub.exp[i],])- colMeans(x[y==sub.exp[i],]))
    }
    if(length(x[y==sub.exp[i],1])==1)
    {
      x[y==sub.exp[i],] <- NA
    }
  }
  return(x)
}


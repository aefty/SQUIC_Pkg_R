# SECTION 4.1: Tests with Synthetic Data
    
# ================================== #
#    Load the necessary packages     #
# ================================== #

library(SQUIC)

if (!require("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}
library(Matrix)

if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

if (!require("MLmetrics", quietly = TRUE)) {
  install.packages("MLmetrics")
}
library(MLmetrics)

if (!require("glasso", quietly = TRUE)) {
  install.packages("glasso")
}
library(glasso)


if (!require("BigQuic", quietly = TRUE)) {
  install.packages("BigQuic")
}
library(BigQuic)

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

if (!require("EQUAL", quietly = TRUE)) {
  devtools::install_github("cescwang85/EQUAL")
}
library(EQUAL)

if (!require("scales", quietly = TRUE)) {
  install.packages("scales")
}
library(scales)

if (!require("latex2exp", quietly = TRUE)) {
  install.packages("latex2exp")
}
library(latex2exp)
# ================================== #


# ================================== #
#    Generate synthetic datasets     #
# ================================== #
generate_data<-function(type,p,n,normalize)
  # type      = "trid" or "rand"
  # p         = number of dimesions
  # n         = number of samples
  # normalize = TRUE or FALSE
{
  set.seed(0);
  
  start_time <- Sys.time()
  
  print(sprintf("# Generating Percision Matrix: type=%s p=%d n=%d",type,p,n));
  
  if(type=="trid") # Tridiagonal matrix for X_star 
  {
    X_star = Matrix::bandSparse(p, p,
                                 (-1):1,
                                 list(rep(-.5, p-1), 
                                      rep(1.25, p), 
                                      rep(-.5, p-1)));
  }
  else if(type=="rand")  # Random matrix for X_star (average of 5 nnz per row) 
  {
    nnz_per_row=5; 
    
    # Make PSD symmetric Random Matrix
    X_star=Matrix::rsparsematrix(nrow=p,ncol=p,nnz=nnz_per_row*p/2,symmetric=TRUE);# we need the divide by 2 (R assuming symmetric)
    x=Matrix::rowSums(abs(X_star))+1;
    D=Matrix::Diagonal(p,x);
    X_star=X_star+D;
    
  }else{
    stop("Unknown matrix type.")
  }
  
  # Generate data
  z    = replicate(n,rnorm(p));
  X_L  = chol(X_star);
  data = matrix(solve(X_L,z),p,n);
  
  finish_time = Sys.time()
  print(sprintf("# Generating Data: time=%f",finish_time-start_time));
  
  if(normalize){
    output <- list(
      "data"   = t(scale(t(data))),
      "X_star" = X_star,
      "var"    = apply(data,1,var)
    );
  }else{
    output <- list(
      "data"   = data, 
      "X_star" = X_star
    );
  }
  
  return(output);
}

# =============================================== #
#  Compare precision matrix estimation packages   #
# =============================================== #
compare <- function(alg,data,lambda,tol,max_iter, X_star, M) 
  # alg      = "SQUIC" or "glasso" or "EQUAL" or "BigQUIC"
  # data     = input data in the form p (dimensions) x n (samples)
  # lambda   = scalar sparsity parameter
  # tol      = termination tolerance
  # max_iter = maximum number of iterations
  # X_star   = true precision matrix
  # M        = graphical bias for the matrix tuning parameter
{
  data_t <- t(data);
  
  verbose = 1;
  
  time_start <- Sys.time()
  
  if(alg=="SQUIC")
  {
    print("#SQUIC")
    # SQUIC
    out	<-SQUIC(
      Y=data,
      lambda=lambda,
      max_iter=max_iter, 
      inv_tol=tol, 
      term_tol=tol, 
      verbose=verbose,
      M=M, 
      X0=NULL, 
      W0=NULL
    );
    
    X	<-out$X;
  }
  else if(alg=="glasso")
  {
    print("#glasso")
    # glasso
    out	<-glasso(
      s=cov(data_t), 
      rho = lambda, 
      nobs=NULL, 
      zero=NULL, 
      thr=tol, 
      maxit=max_iter,  
      approx=FALSE,
      penalize.diagonal=TRUE,w.init=NULL,wi.init=NULL, trace=FALSE);
    X	<- as(out$wi, "sparseMatrix") ;
  }
   else if(alg=="EQUAL")
   {
     print("#EQUAL")
     # EQUAL
     out <- EQUAL(
       X = data_t,
       type = TRUE,
       sdiag = FALSE,
       lambda = lambda,
       lambda.min = sqrt(log(ncol(data_t))/nrow(data_t)),
       nlambda = 1,
       err = tol,
       maxIter = max_iter,
       rho = 1);
     
     X	<- out$Omega[[1]];
   }
  else if(alg=="BigQUIC")
  {
    print("#BigQUIC")
    # BigQUIC
    out	<- BigQuic(
      X = data_t, 
      inputFileName = NULL, 
      outputFileName = NULL, 
      lambda = lambda,
      numthreads = 4, 
      maxit = max_iter, 
      epsilon = tol, 
      k = 0,
      memory_size = 8000,
      verbose = 0, 
      isnormalized = 1, 
      seed = NULL, 
      use_ram = TRUE);
    
    X	<- out$precision_matrices[[1]];
  }
  else
  {
    stop("Alg not found");
  };
  
  time_end	<- Sys.time()
  
  if(!is.null(X_star))
  {
    # Convert matrix to labels
    print("#Computing F1-Score & Accuracy")
    
    X_star_dense_struc = (as.vector(X_star)!=0)*1;
    X_dense_struc = (as.vector(X)!=0)*1;
    
    output <- list(
      "time" = time_end-time_start,
      "X"    = X,
      "fro"  = Matrix::norm(X_star-X,"f")/Matrix::norm(X_star,"f"),	
      "f1"   = MLmetrics::F1_Score(X_star_dense_struc,X_dense_struc,positive="1")			
    );
    
  }else{
    
    # Convert matrix to labels
    output <- list(
      "X"    = X, 
      "time" = time_end-time_start	
    );
  }
  
  return(output);
}


# ================================================================== #
#  Compare the performance of precision matrix estimation packages   #
# ================================================================== #
performance <- function(type,p_set,lambda,n,tol,max_iter)
  # type     = "trid" or "rand"
  # p_set    = the dimensions of the datasets in question
  # lambda   = scalar tuning parameter
  # n        = number of samples
  # tol      = termination tolerance of the algorithms
  # max_iter = maximum number of iterations for all algorithms
{
  l=length(p_set)
  out_squic    =replicate(l, 0);
  out_equal    =replicate(l, 0);	
  out_glasso   =replicate(l, 0);
  out_BigQUIC  =replicate(l, 0);	
  
  for (i in 1:l) {
    
    p=p_set[i];
    
    # Generate the data
    out    = generate_data(type=type, p=p, n=n, normalize = TRUE);
    X_star = out$X_star;
    data   = out$data;
    
    print(sprintf("Benchmark for p=%d started",p));
    
    out = compare(alg="SQUIC", data=data, lambda=lambda, tol=tol, max_iter=max_iter, X_star=NULL, M = NULL);
    out_squic[i] = out$time;
    
    out = compare(alg="EQUAL", data=data, lambda=lambda, tol=tol, max_iter=max_iter, X_star=NULL);
    out_equal[i]=out$time;
    
    out = compare(alg="glasso", data=data, lambda=lambda , tol=tol , max_iter=max_iter , X_star=NULL, M = NULL);
    out_glasso[i] = out$time;
    # 
    out = compare(alg="BigQUIC", data=data, lambda=lambda, tol=tol, max_iter=max_iter, X_star=NULL, M = NULL);
    out_BigQUIC[i] = out$time;
    
  }
  
  output <- list(
    "time_squic"   			= out_squic, 
    "time_equal" 		        = out_equal, 		
    "time_glasso" 			= out_glasso,
    "time_BigQUIC" 			= out_BigQUIC
  )
  
  return(output);
}


# ================================================================== #
#  Evaluate the accuracy in the estimation of precision matrices     #
# ================================================================== #
accuracy <- function(type,lambda_set,p,n,tol,max_iter) 
  # type       = "trid" or "rand"
  # lambda_set = path for the scalar tuning parameter
  # p          = the dimensions of the dataset in question
  # n          = number of samples
  # tol        = termination tolerance of the algorithms
  # max_iter   = maximum number of iterations for all algorithms
{
  
  # Generate data
  out    = generate_data(type=type, p=p, n=n, normalize=TRUE);
  X_star = out$X_star;
  data   = out$data;
  
  l = length(lambda_set)
  out_squic		= replicate(l, 0);
  out_equal		= replicate(l, 0);	
  out_glasso	= replicate(l, 0);
  out_BigQUIC	= replicate(l, 0);
  
  for (i in 1:l) {
    
    lambda = lambda_set[i];
    
    print(sprintf("Benchmark for lambda=%f started",lambda));
    
    out = compare(alg="SQUIC", data=data, lambda=lambda, tol=tol, max_iter=max_iter, X_star=X_star, M = NULL);
    out_squic[i] = out$f1;
    
    out = compare(alg="EQUAL", data=data, lambda=lambda, tol=tol, max_iter=max_iter, X_star=X_star);
    out_equal[i] = out$f1;
    
    out = compare(alg="glasso", data=data, lambda=lambda, tol=tol, max_iter=max_iter, X_star=X_star, M = NULL);
    out_glasso[i] = out$f1;	
    
    out = compare(alg="BigQUIC", data=data, lambda=lambda , tol=tol, max_iter=max_iter, X_star=X_star, M = NULL);
    out_BigQUIC[i]=out$f1;	
  }
  
  output <- list(
    "f1_squic"   			= out_squic, 
    "f1_equal" 		        = out_equal, 		
    "f1_glasso" 			= out_glasso,	
    "f1_BigQUIC" 			= out_BigQUIC	
  )
  
  return(output);
}

# ================================================================== #
#  Evaluate the accuracy in the estimation of precision matrices     #
#  when using a matrix tuning parameter                              #
# ================================================================== #
accuracy_M <- function(type,c,lambda,alpha_set,p,n,tol,max_iter)
  # type      = "trid" or "rand"
  # c         = noise level (additional nonzeros)
  # lambda    = scalar tuning parameter
  # alpha_set = \eta parameter in the estimation of the matrix M
  # p         = the dimensions of the dataset in question
  # n         = number of samples
  # tol       = termination tolerance of the algorithms
  # max_iter  = maximum number of iterations for all algorithms
{
  
  # Generate data
  out     = generate_data( type=type, p=p, n=n, normalize = TRUE);
  X_star  = out$X_star;
  data    = out$data;
  
  # Set the structure of M to be that of X_star .... that would be nice:)
  if(c>0){
    nnz_per_row = c*Matrix::nnzero(X_star)/p;
    M = X_star  + Matrix::rsparsematrix(nrow=p,ncol=p,nnz=nnz_per_row*p/2,symmetric=TRUE);
  }else{
    M = X_star;
  }
  
  l	= length(alpha_set)
  f1_squic	= replicate(l, 0)
  fro_squic	= replicate(l, 0)	
  
  for (i in 1:l) {
    
    alpha = alpha_set[i];
    
    # reset M keeping the structure
    M@x = M@x*0+1
    M   = alpha*M
    
    print(sprintf("Benchmark for alpha=%f started",alpha));
    
    out = compare(alg="SQUIC", data=data, lambda=lambda, tol=tol, max_iter=max_iter, X_star=X_star, M = M)
    f1_squic[i]  = out$f1
    fro_squic[i] = out$fro
  }
  
  output <- list(
    "f1_squic"	= f1_squic,
    "fro_squic"	= fro_squic			
  )
  
  return(output);
}

# =========================== #
# Run the timings experiments #
# =========================== #

# dimensionality of the data
p_set = c(4,16,64,256,1024,4096)

timings_trid = performance(type="trid", p_set=p_set, lambda=0.4, n=100, tol=1e-3, max_iter=100)
timings_rand = performance(type="rand", p_set=p_set, lambda=0.4, n=100, tol=1e-3, max_iter=100)

TRID_ALL_TIMES = cbind(as.data.frame(timings_trid),p_set)
RAND_ALL_TIMES = cbind(as.data.frame(timings_rand),p_set)

# Open the pdf file
pdf("Figure1a.pdf") 
ggplot(TRID_ALL_TIMES, aes(x=p_set)) +
  geom_line(aes(y  = time_squic, colour = "SQUIC"), linetype="solid", alpha=.5) +
  geom_point(aes(y = time_squic, colour = "SQUIC")) +
  geom_line(aes(y  = time_equal, colour = "EQUAL"), linetype="dotdash", alpha=.5) +
  geom_point(aes(y = time_equal, colour = "EQUAL")) +
  geom_line(aes(y  = time_glasso, colour = "glasso"), linetype="twodash", alpha=.5)  +
  geom_point(aes(y = time_glasso, colour = "glasso")) +
  geom_line(aes(y  = time_BigQUIC, colour = "BigQUIC"), linetype="longdash", alpha=.5)  +
  geom_point(aes(y = time_BigQUIC, colour = "BigQUIC")) +
  scale_x_continuous(breaks=p_set, trans = "log10") +
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  ylab(TeX("Runtime (sec)")) +
  xlab(TeX("Dimensions $(p)$")) +
  ggtitle(TeX("Figure 1a: Timings for the tridiagonal dataset")) +
  theme_bw() +
  theme(legend.title=element_blank())
# Close the pdf file
dev.off()


# Open the pdf file
pdf("Figure1b.pdf") 
ggplot(RAND_ALL_TIMES, aes(x=p_set)) +
  geom_line(aes(y  = time_squic, colour = "SQUIC"), linetype="solid", alpha=.5) +
  geom_point(aes(y = time_squic, colour = "SQUIC")) +
  geom_line(aes(y  = time_equal, colour = "EQUAL"), linetype="dotdash", alpha=.5) +
  geom_point(aes(y = time_equal, colour = "EQUAL")) +
  geom_line(aes(y  = time_glasso, colour = "glasso"), linetype="twodash", alpha=.5)  +
  geom_point(aes(y = time_glasso, colour = "glasso")) +
  geom_line(aes(y  = time_BigQUIC, colour = "BigQUIC"), linetype="longdash", alpha=.5)  +
  geom_point(aes(y = time_BigQUIC, colour = "BigQUIC")) +
  scale_x_continuous(breaks=p_set, trans = "log10") +
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  ylab(TeX("Runtime (sec)")) +
  xlab(TeX("Dimensions $(p)$")) +
  ggtitle(TeX("Figure 1b: Timings for the random dataset")) +
  theme_bw() +
  theme(legend.title=element_blank())
# Close the pdf file
dev.off()

# ============================================================  #
# Run the accuracy experiments with a scalar tuning parameter   #
# ============================================================  #

# lambda path 
lambda_set = c(.2,.25,.3,.35,.4,.45,.5,.55,.6)

accuracy_trid = accuracy(type="trid",lambda_set=lambda_set,p=1024,n=100,tol=1e-4,max_iter=100) 
accuracy_rand = accuracy(type="rand",lambda_set=lambda_set,p=1024,n=100,tol=1e-4,max_iter=100) 

TRID_ALL_FSCORE = cbind(as.data.frame(accuracy_trid),lambda_set)
RAND_ALL_FSCORE = cbind(as.data.frame(accuracy_rand),lambda_set)

# Open the pdf file
pdf("Figure2a.pdf") 
ggplot(TRID_ALL_FSCORE, aes(x=lambda_set)) +
  geom_line(aes(y  = f1_squic, colour = "SQUIC"), linetype="solid", alpha=.5) +
  geom_point(aes(y = f1_squic, colour = "SQUIC")) +
  geom_line(aes(y  = f1_equal, colour = "EQUAL"), linetype="dotdash", alpha=.5) +
  geom_point(aes(y = f1_equal, colour = "EQUAL")) +
  geom_line(aes(y  = f1_glasso, colour = "glasso"), linetype="twodash", alpha=.5)  +
  geom_point(aes(y = f1_glasso, colour = "glasso")) +
  geom_line(aes(y  = f1_BigQUIC, colour = "BigQUIC"), linetype="longdash", alpha=.5)  +
  geom_point(aes(y = f1_BigQUIC, colour = "BigQUIC")) +
  ylab(TeX("F1-Score")) +
  xlab(TeX("Scalar Sparsity Parameter $(\\lambda)$")) +
  ggtitle(TeX("Figure 2a: Accuracy for the tridiagonal dataset")) +
  theme_bw() +
  theme(legend.title=element_blank())
# Close the pdf file
dev.off()

# Open the pdf file
pdf("Figure2b.pdf") 
ggplot(RAND_ALL_FSCORE, aes(x=lambda_set)) +
  geom_line(aes(y  = f1_squic, colour = "SQUIC"), linetype="solid", alpha=.5) +
  geom_point(aes(y = f1_squic, colour = "SQUIC")) +
  geom_line(aes(y  = f1_equal, colour = "EQUAL"), linetype="dotdash", alpha=.5) +
  geom_point(aes(y = f1_equal, colour = "EQUAL")) +
  geom_line(aes(y  = f1_glasso, colour = "glasso"), linetype="twodash", alpha=.5)  +
  geom_point(aes(y = f1_glasso, colour = "glasso")) +
  geom_line(aes(y  = f1_BigQUIC, colour = "BigQUIC"), linetype="longdash", alpha=.5)  +
  geom_point(aes(y = f1_BigQUIC, colour = "BigQUIC")) +
  ylab(TeX("F1-Score")) +
  xlab(TeX("Scalar Sparsity Parameter $(\\lambda)$")) +
  ggtitle(TeX("Figure 2b: Accuracy for the random dataset")) +
  theme_bw() +
  theme(legend.title=element_blank())
# Close the pdf file
dev.off()

# ============================================================  #
# Run the accuracy experiments with a matrix tuning parameter   #
# ============================================================  #

# level of noise
c1 = 0; c2 = 2; c3 = 10;
# scalar sparsity parameter
lambda =.95
# \eta parameter for the matrix M
alpha_set = c(.0001,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)

# testing the tridiagonal dataset
accuracy_M_trid_1 = accuracy_M(type = "trid", c = c1, lambda = lambda,
                              alpha_set = alpha_set, p = 1024, n = 100, tol=1e-3,max_iter=100)

accuracy_M_trid_2 = accuracy_M(type = "trid", c = c2, lambda = lambda,
                                alpha_set = alpha_set, p = 1024, n = 100, tol=1e-3,max_iter=100)

accuracy_M_trid_3 = accuracy_M(type = "trid", c = c3, lambda = lambda,
                                alpha_set = alpha_set, p = 1024, n = 100, tol=1e-3,max_iter=100)

FSCORE_M_TRID = cbind(as.data.frame(accuracy_M_trid_1$f1_squic), as.data.frame(accuracy_M_trid_2$f1_squic),
                      as.data.frame(accuracy_M_trid_3$f1_squic), alpha_set)


# Open the pdf file
pdf("Figure3a.pdf") 
ggplot(FSCORE_M_TRID, aes(x=alpha_set)) +
  geom_line(aes(y  = accuracy_M_trid_1$f1_squic, colour = "c = 0"), linetype="solid", alpha=.5) +
  geom_point(aes(y = accuracy_M_trid_1$f1_squic, colour = "c = 0")) +
  geom_line(aes(y  = accuracy_M_trid_2$f1_squic, colour = "c = 2"), linetype="twodash", alpha=.5)  +
  geom_point(aes(y = accuracy_M_trid_2$f1_squic, colour = "c = 2")) +
  geom_line(aes(y  = accuracy_M_trid_3$f1_squic, colour = "c = 10"), linetype="longdash", alpha=.5)  +
  geom_point(aes(y = accuracy_M_trid_3$f1_squic, colour = "c = 10")) +
  ylab(TeX("F1-Score")) +
  xlab(TeX("Bias Parameter $(\\eta)$")) +
  ggtitle(TeX("Figure 3a: Accuracy for the tridiagonal dataset")) +
  theme_bw() +
  labs(color='Noise Level')
# Close the pdf file
dev.off()


# testing the random dataset
accuracy_M_rand_1 = accuracy_M(type = "rand", c = c1, lambda = lambda,
                                alpha_set = alpha_set, p = 1024, n = 100, tol=1e-3,max_iter=100)

accuracy_M_rand_2 = accuracy_M(type = "rand", c = c2, lambda = lambda,
                                alpha_set = alpha_set, p = 1024, n = 100, tol=1e-3,max_iter=100)

accuracy_M_rand_3 = accuracy_M(type = "rand", c = c3, lambda = lambda,
                                alpha_set = alpha_set, p = 1024, n = 100, tol=1e-3,max_iter=100)

FSCORE_M_RAND = cbind(as.data.frame(accuracy_M_rand_1$f1_squic), as.data.frame(accuracy_M_rand_2$f1_squic),
                      as.data.frame(accuracy_M_rand_3$f1_squic), alpha_set)


# Open the pdf file
pdf("Figure3b.pdf") 
ggplot(FSCORE_M_RAND, aes(x=alpha_set)) +
  geom_line(aes(y  = accuracy_M_rand_1$f1_squic, colour = "c = 0"), linetype="solid", alpha=.5) +
  geom_point(aes(y = accuracy_M_rand_1$f1_squic, colour = "c = 0")) +
  geom_line(aes(y  = accuracy_M_rand_2$f1_squic, colour = "c = 2"), linetype="twodash", alpha=.5)  +
  geom_point(aes(y = accuracy_M_rand_2$f1_squic, colour = "c = 2")) +
  geom_line(aes(y  = accuracy_M_rand_3$f1_squic, colour = "c = 10"), linetype="longdash", alpha=.5)  +
  geom_point(aes(y = accuracy_M_rand_3$f1_squic, colour = "c = 10")) +
  ylab(TeX("F1-Score")) +
  xlab(TeX("Bias Parameter $(\\eta)$")) +
  ggtitle(TeX("Figure 3b: Accuracy for the random dataset")) +
  theme_bw() +
  labs(color='Noise Level')
# Close the pdf file
dev.off()

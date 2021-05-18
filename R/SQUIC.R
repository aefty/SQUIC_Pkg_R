usethis::use_package("Matrix") 

# Here we set the envioerment variable KMP_DUPLICATE_LIB_OK=TRUE
# This is needed for Mac due to potential conflict with other OMP versions
#.onLoad <- function(libname, pkgname){
#	print("## SQUIC init start ##");
#	print("1) Setting envirnment variable KMP_DUPLICATE_LIB_OK=TRUE")
#	Sys.setenv(KMP_DUPLICATE_LIB_OK =TRUE)
#	print("## SQUIC init finished ##");
#}

# Main function
SQUIC <- function(Y, lambda, max_iter, drop_tol, term_tol,verbose=1, M=NULL, X0=NULL, W0=NULL) {
  
  
  verbose = min(verbose,1);
  
  Y = as.matrix(Y);
  p=nrow(Y);
  n=ncol(Y);

  if(p<3){
	  stop('#SQUIC: number of random variables (p) must larger than 2');
  }
  
  if(n<2){
    stop('#SQUIC: number of samples (n) must be larger than 1 .');
  }
  
  if(lambda<=0){
	  stop('#SQUIC: lambda must great than zero.');
  }
  if(max_iter<0){
	  stop('#SQUIC: max_iter cannot be negative.');
  }
  
  if(drop_tol<=0){
	  stop('#SQUIC: drop_tol must great than zero.');
  }
  
  if(term_tol<=0){
	  stop('#SQUIC: term_tol must great than zero.');
  } 

 
  if(is.null(M)){
	  # Make empty sparse matrix of type dgCMatrix.
	  M	= as(Matrix::sparseMatrix(dims = c(p,p), i={}, j={}),"dgCMatrix");
  }else{
    # Force symmetrix
    M<-abs(M);
  	M<-Matrix::forceSymmetric(M,uplo="L"); 
  }

  if(is.null(X0) || is.null(W0)){
	  # Make empty sparse matrix of type dgCMatrix.
	  X0 = as(Matrix::sparseMatrix(dims = c(p,p), i=c(1:p), j=c(1:p) , x=rep(1, p)),"dgCMatrix");
	  W0 = as(Matrix::sparseMatrix(dims = c(p,p), i=c(1:p), j=c(1:p) , x=rep(1, p)),"dgCMatrix");
  }else{
    # Force symmetrix
    X0<-Matrix::forceSymmetric(X0,uplo="L"); 
    W0<-Matrix::forceSymmetric(W0,uplo="L"); 
  }
  
  #Use block variant of SQUIC
  mode=0; 

  output = SQUIC::SQUIC_R(
   	Y , 
   	lambda , 
   	max_iter , drop_tol , term_tol , 
   	verbose , mode , 
   	M , X0 , W0);

  return(output);
}

# Sparse Sample covariance matrix
SQUIC_S<-function(Y, lambda,verbose=1, M=NULL){

	# Get sample covarinace matrix by running SQUIC with max_iter=0;
	output = SQUIC::SQUIC( 
		Y=Y , 
		lambda=lambda ,
		max_iter=0 , drop_tol=1e-6 , term_tol=1e-6 , 
		verbose=verbose , M=M , X0=NULL , W0=NULL);

	return(output);
}

# Sparse Sample covariance matrix
SQUIC_CV<-function(Y, lambda_set, M=NULL, X0=NULL, W0=NULL){
  
  l<-length(lambda_set)
  klc		<-replicate(l, 0);
  tttt		<-replicate(l, 0);
  
  tol = (1e-3);
  
  for (i in 1:l) {
    
    lambda= lambda_set[i];
    
    out<-SQUIC::SQUIC(Y=Y, lambda=lambda, max_iter=5, drop_tol=tol/2, term_tol=tol,verbose=1, M=M, X0=X0, W0=X0);
    
    ll = (n/2)*(out$info_logdetX-out$info_trSX)
    
    
    
    R = t(Y)%*%out$X%*%(Y);
    
    temp=0; # sum(diag(R)^2)
    for (j in 1:n) {
      temp = temp+R[j,j]^2;
    }
    

    
    klc[i] = -(1/n)*ll + 1/(2*n*(n-1) ) * (temp - n*sum(R^2)/n^2);
    
    
    
    
    
    
    tttt[i]=n*sum(R^2)/n^2
    
  }
  output <- list(
    "scores" = klc,
    "lambda_opt"    = lambda_set[which.min(klc)]
  );
  
  
  print(tttt)
  
  
  return(output);
}




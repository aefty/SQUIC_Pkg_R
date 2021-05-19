usethis::use_package("Matrix") 

# Main function
SQUIC <- function(Y, lambda, max_iter=100, drop_tol=1e-4, term_tol=1e-3,verbose=1, M=NULL, X0=NULL, W0=NULL) {
  
  verbose <- min(verbose,1);
  
  Y <- as.matrix(Y);
  p <- nrow(Y);
  n <- ncol(Y);

  if(p<3){
	  stop('#SQUIC: number of random variables (p) must larger than 2');
  }
  if(n<2){
    stop('#SQUIC: number of samples (n) must be larger than 1 .');
  }
  if(lambda<=0){
	  stop('#SQUIC: lambda must be great than zero.');
  }
  if(max_iter<0){
	  stop('#SQUIC: max_iter cannot be negative.');
  }
  if(drop_tol<=0){
	  stop('#SQUIC: drop_tol must be great than zero.');
  }
  if(term_tol<=0){
	  stop('#SQUIC: term_tol must be great than zero.');
  } 

  if(is.null(M)){
	  # Make empty sparse matrix of type dgCMatrix.
	  M	<- as(Matrix::sparseMatrix(dims = c(p,p), i={}, j={}),"dgCMatrix");
  }else{

    # Make all postive, drop all zeros and force symmetrix
    M <- abs(M);
    M <- drop0(M);
  	M <- Matrix::forceSymmetric(M,uplo="L"); 
  }

  if(is.null(X0) || is.null(W0)){
	  # Make empty sparse matrix of type dgCMatrix.
	  X0 <- as(Matrix::sparseMatrix(dims = c(p,p), i=c(1:p), j=c(1:p) , x=rep(1, p)),"dgCMatrix");
	  W0 <- as(Matrix::sparseMatrix(dims = c(p,p), i=c(1:p), j=c(1:p) , x=rep(1, p)),"dgCMatrix");
  }else{
    # Force symmetrix
    X0 <- Matrix::forceSymmetric(X0,uplo="L"); 
    W0 <- Matrix::forceSymmetric(W0,uplo="L"); 
  }
  
  #Use block variant of SQUIC
  mode <- 0; 

  output <- SQUIC::SQUIC_R(
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
	output <- SQUIC::SQUIC( 
		Y        =Y, 
		lambda   =lambda,
		max_iter =0,  
		verbose  =verbose, 
    M        =M, 
    X0       =NULL,
    W0       =NULL);

	return(output);
}


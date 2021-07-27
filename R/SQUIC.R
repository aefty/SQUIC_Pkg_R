usethis::use_package("Matrix") 

# Main function
SQUIC <- function(Y, lambda, max_iter=100, tol=1e-3,verbose=1, M=NULL, X0=NULL, W0=NULL) {
  
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
  if(tol<=0){
	  stop('#SQUIC: tol must be great than zero.');
  }

  if(is.null(M)){
	  # Make empty sparse matrix of type dgCMatrix.
	  M = as(Matrix::sparseMatrix(dims = c(p,p), i={}, j={}),"dgCMatrix");
  }else{

    if(nrow(M)!= p || ncol(M)!=p){
      stop('#SQUIC: M must be square matrix with size pxp.');
    } 

    # Make all postive, drop all zeros and force symmetric
    M = Matrix::drop0(M);
    M <- abs(M);
  	M <- Matrix::forceSymmetric(M,uplo="L"); 
  }

  if(is.null(X0) || is.null(W0)){
	  # Make identity sparse matrix of type dgCMatrix.
	  X0 <- as(Matrix::sparseMatrix(dims = c(p,p), i=c(1:p), j=c(1:p) , x=rep(1, p)),"dgCMatrix");
	  W0 <- as(Matrix::sparseMatrix(dims = c(p,p), i=c(1:p), j=c(1:p) , x=rep(1, p)),"dgCMatrix");
  }else{

    if(nrow(X0)!= p || ncol(X0)!=p){
      stop('#SQUIC: X0 must be square matrix with size pxp.');
    } 

    if(nrow(W0)!= p || ncol(W0)!=p){
      stop('#SQUIC: W0 must be square matrix with size pxp.');
    } 

    # Force Symmetric
    X0 <- Matrix::forceSymmetric(X0,uplo="L"); 
    W0 <- Matrix::forceSymmetric(W0,uplo="L"); 
  }
  
  #Use block variant of SQUIC
  mode <- 0; 
  
  # Hard code both tollerance to be the same
  term_tol = tol;
  inv_tol  = tol;

  output <- SQUIC::SQUIC_R(
   	Y , 
   	lambda , 
   	max_iter , inv_tol , term_tol , 
   	verbose , mode , 
   	M , X0 , W0);

  return(output);
}

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
SQUIC <- function(Y, lambda, max_iter, drop_tol, term_tol,verbose=1, mode=0, M=NULL, X0=NULL, W0=NULL) {
  
  p=nrow(Y);

  if(p<3){
	  stop('#SQUIC: size of matrix must be larger than 2.');
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

  if(mode < 0 || mode > 9){
	  stop('#SQUIC: mode must be in [0,9].');
  } 

  if(is.null(M)){
	  # Make empty sparse matrix of type dgCMatrix.
	  M	= as(Matrix::sparseMatrix(dims = c(p,p), i={}, j={}),"dgCMatrix");
  }else{
    if(!isSymmetric(M)){
      stop('#SQUIC: M must be symmetric.');
    }
  }

  if(is.null(X0) || is.null(W0)){
	  # Make empty sparse matrix of type dgCMatrix.
	  X0 = as(Matrix::sparseMatrix(dims = c(p,p), i=c(1:p), j=c(1:p) , x=rep(1, p)),"dgCMatrix");
	  W0 = as(Matrix::sparseMatrix(dims = c(p,p), i=c(1:p), j=c(1:p) , x=rep(1, p)),"dgCMatrix");
  }else{
     if(!isSymmetric(X0) || isSymmetric(W0)){
      stop('#SQUIC: X0 and W0 must be symmetric.')
    }
  }

   output = SQUIC::SQUIC_R(
   	Y , 
   	lambda , 
   	max_iter , drop_tol , term_tol , 
   	verbose , mode , 
   	M , X0 , W0);

   return(output);
}

# Sparse Sample covariance matrix
SQUIC_S<-function(Y, lambda, M=NULL,verbose=1){

	# Get sample covarinace matrix by running SQUIC with max_iter=0;
	output = SQUIC::SQUIC( 
		Y=Y , 
		lambda=lambda ,
		max_iter=0 , drop_tol=1e-6 , term_tol=1e-6 , 
		verbose=verbose , mode=0,  
		M=M , X0=NULL , W0=NULL);

	return(output);
}


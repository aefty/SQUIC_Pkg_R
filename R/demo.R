usethis::use_package("MLmetrics") 

DEMO.generate_data<-function(type="trid",p=4^5,n=100,normalize=TRUE)
{
	set.seed(1);

    start_time <- Sys.time()
	
    print(sprintf("# Generating Percision Matrix: type=%s p=%d n=%d",type,p,n));

	if(type=="eye") # Idendity Matrix for iC_star
	{
		iC_star <- Matrix::Diagonal(p);
	}
	else if(type=="trid") # Tridiagiaonl matrix for iC_star 
	{
		iC_star <- Matrix::bandSparse(p, p,
				(-1):1,
				list(rep(-.5, p-1), 
					rep(1.25, p), 
					rep(-.5, p-1)));
	}
		else if(type=="rand")  # Random matrix for iC_star (averag of 5 nnz per row) 
	{
		nnz_per_row=5; 

		# Make PSD symmetric Random Matrix
		iC_star <-Matrix::rsparsematrix(p,p,NULL,nnz_per_row*p/2,symmetric=TRUE);# we need the divide by 2 (R assuming symmetric)
		x=Matrix::rowSums(abs(iC_star))+1;
		D=Matrix::Diagonal(p,x);
		iC_star<- iC_star+D;

	}else{
		stop("Unknown matrix type.")
	}

	# Generate data
	z    <- replicate(n,rnorm(p));
	iC_L <- chol(iC_star);
	data <- matrix(solve(iC_L,z),p,n);

	finish_time <- Sys.time()
	print(sprintf("# Generating Data: time=%f",finish_time-start_time));


	if(normalize){
        output <- list(
            "data"   = t(scale(t(data))),
            "X_star" = iC_star,
            "var"    = apply(data,1,var)
        );
	}else{
		output <- list(
			"data"   = data, 
			"X_star" = iC_star
		);
	}

	return(output);
}


DEMO.performance <- function(type,p_set=c(4,16,64,256,1024,4096),lambda=0.4,n=100,tol=1e-4,max_iter=10) 
{

	l<-length(p_set)
	out_squic		<-replicate(l, 0);
	out_equal		<-replicate(l, 0);	
	out_glasso		<-replicate(l, 0);

	for (i in 1:l) {

        p=p_set[i];

        # Generate data
	    out<-SQUIC::DEMO.generate_data( type=type , p=p , n=n );
	    X_star<-out$X_star;
	    data<-out$data;

		print(sprintf("Benchmark for p=%d started",p));

		out<-SQUIC::DEMO.compare(alg="SQUIC"   , data=data , lambda=lambda , tol=tol , max_iter=max_iter , X_star=NULL);
		out_squic[i]<-out$time;

		out<-SQUIC::DEMO.compare(alg="EQUAL"   , data=data , lambda=lambda , tol=tol , max_iter=max_iter , X_star=NULL);
		out_equal[i]<-out$time;

		out<-SQUIC::DEMO.compare(alg="glasso"    ,  data=data , lambda=lambda , tol=tol , max_iter=max_iter , X_star=NULL);
		out_glasso[i]<-out$time;			
	}

	output <- list(
		"time_squic"   			= out_squic, 
		"time_equal" 			= out_equal, 		
		"time_glasso" 			= out_glasso 			
		)

	return(output);
}


DEMO.CV <- function(type,lambda_set=c(.2,.25,.3,.35,.4,.45,.5,.55,.6),p=1024,n=100,tol=1e-4,max_iter=10) 
{
  
  # Generate data
  out<-SQUIC::DEMO.generate_data( type=type , p=p , n=n );
  X_star<-out$X_star;
  data<-out$data;
  
  print(sprintf("Benchmark for lambda=%f started",lambda));
    
  out<-SQUIC::SQUIC_CV(Y=data,lambda_set = lambda_set);


  return(out);
}


DEMO.GR<- function(type,lambda_set=c(.2,.25,.3,.35,.4,.45,.5,.55,.6),p=1024,n=100,tol=1e-4,max_iter=10) 
{
  
  
}


DEMO.accuracy <- function(type,lambda_set=c(.2,.25,.3,.35,.4,.45,.5,.55,.6),p=1024,n=100,tol=1e-4,max_iter=10) 
{

    # Generate data
	out<-SQUIC::DEMO.generate_data( type=type , p=p , n=n );
	X_star<-out$X_star;
	data<-out$data;

	l<-length(lambda_set)
	out_squic		<-replicate(l, 0);
	out_equal		<-replicate(l, 0);	
	out_glasso		<-replicate(l, 0);

	for (i in 1:l) {

        lambda=lambda_set[i];

		print(sprintf("Benchmark for lambda=%f started",lambda));

		out<-SQUIC::DEMO.compare(alg="SQUIC"   , data=data , lambda=lambda , tol=tol , max_iter=max_iter , X_star=X_star);
		out_squic[i]<-out$f1;

		out<-SQUIC::DEMO.compare(alg="EQUAL"   , data=data , lambda=lambda , tol=tol , max_iter=max_iter , X_star=X_star);
		out_equal[i]<-out$f1;

		out<-SQUIC::DEMO.compare(alg="glasso"  ,  data=data , lambda=lambda , tol=tol , max_iter=max_iter , X_star=X_star);
		out_glasso[i]<-out$f1;			
	}

	output <- list(
		"f1_squic"   			= out_squic, 
		"f1_equal" 			    = out_equal, 		
		"f1_glasso" 			= out_glasso 			
		)

	return(output);
}



DEMO.accuracy_M <- function(type,corrupt_nnz_fac=0,lambda=.5,alpha_set=c(.0,.1,.2,.5,.6,.7,.8,.9,1),p=1024,n=100,tol=1e-4,max_iter=10) 
{

  # Generate data
	out<-SQUIC::DEMO.generate_data( type=type , p=p , n=n );
	X_star<-out$X_star;
	data<-out$data;

	# Set the structure of M to be that of X_star .... that would be nice:)
	if(corrupt_nnz_fac>0){
		nnz_per_row = corrupt_nnz_fac*Matrix::nnzero(X_star)/p;
		M = X_star  + Matrix::rsparsematrix(p,p,NULL,nnz_per_row*p/2,symmetric=TRUE);
	}else{
		M = X_star;
	}

	l           <-length(alpha_set)
	f1_squic		<-replicate(l, 0);
	fro_squic		<-replicate(l, 0);	

	for (i in 1:l) {
	  
	  alpha=alpha_set[i];
	  
	  # reset M keeping the structure
	  M@x=M@x*0+1;
	  M = alpha*lambda*M;
  
		print(sprintf("Benchmark for alpha=%f started",alpha));

		out<-SQUIC::DEMO.compare(alg="SQUIC" , data=data , lambda=lambda , tol=tol , max_iter=max_iter , X_star=X_star, M = M );
		f1_squic[i]<-out$f1;
		fro_squic[i]<-out$fro;
	}

	output <- list(
		"f1_squic"	= f1_squic,
		"fro_squic"	= fro_squic			
	)

	return(output);
}



DEMO.compare <- function(alg,data,lambda=0.5,tol=1e-4,max_iter=10, X_star= NULL, M=NULL) 
{
	data_t <- Matrix::t(data);

	verbose = 1;

	time_start <- Sys.time()

	if(alg=="SQUIC")
	{
		print("#SQUIC")
		# SQUIC
		out	<-SQUIC::SQUIC(
			Y=data,
			lambda=lambda,
			max_iter=max_iter, 
			drop_tol=tol/2, 
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
		out	<-glasso::glasso(
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
		out	<-EQUAL::EQUAL(
			X=data_t,
			type=TRUE,
			sdiag=FALSE,
			lambda=lambda,
			lambda.min=sqrt(log(ncol(data_t))/nrow(data_t)),
			nlambda=1,
			err=tol,
			maxIter = max_iter,
			rho=1);

		X	<-out$Omega[[1]];
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

\name{SQUIC}
\alias{SQUIC}
\docType{package}

\title{
	\packageTitle{SQUIC}
}

\description{
SQUIC is a second-order, L1-regularized maximum likelihood method for performant large-scale sparse precision matrix estimation.
}

\usage{
 SQUIC(Y, lambda, max_iter = 100, tol = 1e-3, verbose = 1, M = NULL,  X0=NULL, W0 = NULL)
 }
 \arguments{
 \item{Y}{input data in the form p (dimensions) x n (samples)}
 \item{lambda}{(Non-negative) Scalar tuning parameter controlling the sparsity of the precision matrix. lambda=0 means no regularization.}
 \item{max_iter}{Maximum number of Newton iterations of outer loop. Default: 100}
 \item{tol}{Tolerance for convergence. Iterations stop when ... Default: 1e-3}
 \item{verbose}{Level of printing output. Default: 1}
 \item{X0}{Optional starting values for the estimated precision matrix (p by p). Default: NULL}
 \item{wi.init}{Optional starting values for the estimated inverse of the precision matrix (p by p).
 Default: NULL}
 }


\details{
 Estimates a sparse precision matrix using the L1-regularized maximum likelihood (ML) 
 approach, commonly referred to as the "graphical lasso" method. 
 }
 
 \value{
 A list with components
 \item{X}{Estimated precision matrix}
 \item{W}{Estimated inverse of the precision matrix}
 \item{info_time_total}{Total time spent}
 \item{info_time_sample_cov}{Time spent for the estimation of the sample covariance matrix}
 \item{info_time_optimize}{Time spent in the ...}
 \item{info_time_factor}{Time spent in the Cholesky factorization}
 \item{info_time_approximate_inv}{Time spent in the approximate inversion of X}
 \item{info_time_coordinate_upd}{Time spend in the block coordinate descent update}
 \item{info_objective}{Value of the objective at each Newton iteration}
 \item{info_logdetX}{Value of the log determinant of X}
 \item{info_trSX}{Value of trace(SX)}
 }

\references{
Bollhöfer, M., Eftekhari, A., Scheidegger, S. and Schenk, O., 2019. Large-scale sparse inverse covariance   matrix estimation. SIAM Journal on Scientific Computing, 41(1), pp.A380-A401.

Eftekhari, A., Bollhöfer, M. and Schenk, O., 2018, November. Distributed memory sparse inverse covariance      matrix estimation on high-performance computing architectures. In SC18: International Conference for High   Performance Computing, Networking, Storage and Analysis (pp. 253-264). IEEE.

Eftekhari A, Pasadakis D, Bollhöfer M, Scheidegger S, Schenk O (2021). Block-Enhanced PrecisionMatrix          Estimation for Large-Scale Datasets. Journal of Computational Science, p. 101389.

\examples{

 library(SQUIC)
 library(Matrix)

 p=10
 n=130
 lambda=.5
 max_iter=10
 tol=1e-3

 # generate a tridiagonal matrix
 iC_star = Matrix::bandSparse(p, p, (-1):1, list(rep(-.5, p-1), rep(1.25, p), rep(-.5, p-1)));

 # generate the data
 z    = replicate(n,rnorm(p));
 iC_L = chol(iC_star);
 data = matrix(solve(iC_L,z),p,n);

# Run SQUIC
 out<-SQUIC(data,lambda)
 }


 \keyword{covariance matrices}
 \keyword{multivariate}
 \keyword{models}
 \keyword{graphs}

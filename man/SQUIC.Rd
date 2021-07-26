\name{SQUIC}
\alias{SQUIC}

\title{
	\packageTitle{SQUIC}
}

\description{
SQUIC is a second-order, L1-regularized maximum likelihood method for performant large-scale sparse precision matrix estimation.
}

\usage{
 SQUIC(Y, lambda, max_iter = 100, tol = 1e-3, verbose = 1, M = NULL,  X0 = NULL, W0 = NULL)
 }
 \arguments{
 \item{Y}{Input data in the form p (dimensions) by n (samples).}
 \item{lambda}{(Non-zero positive paramater) Scalar tuning parameter controlling the sparsity of the precision matrix.}
 \item{max_iter}{Maximum number of Newton iterations of outer loop. Default: 100.}
 \item{tol}{Tolerance for convergence and approximate inversion. Default: 1e-3.}
 \item{verbose}{Level of printing output (0 or 1). Default: 1}
 \item{M}{The matrix encoding the sparsity pattern of the matrix tuning parameter, i.e., \Lambda (p by p). Default: NULL.}
 \item{X0}{Initial guess of precision matrix (p by p). Default: NULL.}
 \item{W0}{Initial guess of the inverse of the precision matrix (p by p). Default: NULL.}
 }


\details{
SQUIC is a high-performance precision matrix estimation package intended for large-scale applications. The algorithm utilizes second-order information in conjunction with a highly optimized, supernodal Cholesky factorization routine (CHOLMOD), a block-oriented computational approach for both the coordinate descent and the approximate matrix inversion, and finally, a sparse representation of the sample covariance matrix. The library is implemented in C++ with shared-memory parallelization using OpenMP.
For further details, please see the listed references.
}
 
 \value{
 A list with components
 \item{X}{Estimated precision matrix (p by p).}
 \item{W}{Estimated inverse of the precision matrix (p by p).}
 \item{info_time_total}{Total runtime.}
 \item{info_time_sample_cov}{Runtime of the sample covariance matrix.}
 \item{info_time_optimize}{Runtime of the Newton steps.}
 \item{info_time_factor}{Runtime of the Cholesky factorization.}
 \item{info_time_approximate_inv}{Runtime of the approximate matrix inversion.}
 \item{info_time_coordinate_upd}{Runtime of the coordinate descent update.}
 \item{info_objective}{Value of the objective at each Newton iteration.}
 \item{info_logdetX}{Value of the log determinant of the precision matrix.}
 \item{info_trSX}{Value of the trace of the sample covariance matrix times the precision matrix.}
 }

\references{
Bollhöfer, M., Eftekhari, A., Scheidegger, S. and Schenk, O., 2019. Large-scale sparse inverse covariance matrix estimation. SIAM Journal on Scientific Computing, 41(1), pp.A380-A401.

Eftekhari, A., Bollhöfer, M. and Schenk, O., 2018, November. Distributed memory sparse inverse covariance matrix estimation on high-performance computing architectures. In SC18: International Conference for High Performance Computing, Networking, Storage and Analysis (pp. 253-264). IEEE.

Eftekhari, A., Pasadakis, D., Bollhöfer, M., Scheidegger, S. and Schenk, O., 2021. Block-Enhanced Precision Matrix          Estimation for Large-Scale Datasets. Journal of Computational Science, p. 101389.
}


\examples{

 library(SQUIC)

 p = 1024
 n = 100
 lambda = .4
 max_iter = 100
 tol = 1e-3

 # generate a tridiagonal matrix
 iC_star = Matrix::bandSparse(p, p, (-1):1, list(rep(-.5, p-1), rep(1.25, p), rep(-.5, p-1)));

 # generate the data
 z    = replicate(n,rnorm(p));
 iC_L = chol(iC_star);
 data = matrix(solve(iC_L,z),p,n);

 # Run SQUIC
 out <- SQUIC(data,lambda)
 }

 \keyword{precision matrix estimation}
 \keyword{multivariate analysis}
 \keyword{graphical lasso}

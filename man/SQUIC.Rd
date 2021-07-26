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
 glasso(s, rho, nobs=NULL, zero=NULL, thr=1.0e-4, maxit=1e4,  approx=FALSE,
 penalize.diagonal=TRUE, start=c("cold","warm"),
 w.init=NULL,wi.init=NULL, trace=FALSE)
 }
 \arguments{
 \item{s}{Covariance matrix:p by p matrix (symmetric)}
 \item{rho}{(Non-negative) regularization parameter for lasso. rho=0 means no regularization. Can be a scalar    (usual) or a symmetric p by p matrix, or a vector of length p. In the latter case, the penalty matrix has jkth  element sqrt(rho[j]*rho[k]).}
 \item{nobs}{Number of observations used in computation of the  covariance matrix s. This quantity is need to    compute the value of  log-likelihood.
 If not specified, loglik will be returned as NA.}
 \item{zero}{(Optional) indices of entries of inverse covariance to be constrained to be zero. The input should  be a matrix with two columns, each row indicating
 the indices of elements to be constrained to be zero. The solution must be symmetric, so you need only          specify one of (j,k) and (k,j). An entry in the zero matrix
 overrides any entry in the rho matrix for a given element.}
 \item{thr}{Threshold for convergence. Default value is 1e-4.  Iterations stop when average absolute parameter   change is less than thr * ave(abs(offdiag(s)))}
 \item{maxit}{Maximum number of iterations of outer loop. Default 10,000}
 \item{approx}{Approximation flag: if true, computes Meinhausen-Buhlmann(2006)
 approximation}
 \item{penalize.diagonal}{Should diagonal of inverse covariance be penalized?
 Dafault TRUE.}
 \item{start}{Type of start. Cold start is default. Using Warm start, can provide starting values for w and wi}
 \item{w.init}{Optional starting values for estimated covariance matrix (p by p).
 Only needed when start="warm" is specified}
 \item{wi.init}{Optional starting values for estimated inverse covariance matrix (p by p)
 Only needed when start="warm" is specified}
 \item{trace}{Flag for printing out information as iterations proceed.
 Default FALSE}
 }


\details{
 Estimates a sparse inverse covariance matrix using a lasso (L1) penalty,
 using the approach of Friedman, Hastie and Tibshirani (2007).
 The Meinhausen-Buhlmann (2006) approximation is also implemented.
 The algorithm can also be used to estimate a graph with missing edges,
 by specifying which edges to omit in the zero argument, and setting rho=0.
 Or both fixed zeroes for some elements and regularization on the other elements
 can be specified.

 This version 1.7 uses a block diagonal screening rule to speed up
 computations considerably. Details are given in the paper "New insights
 and fast computations for the graphical lasso" by Daniela Witten, Jerry
 Friedman, and Noah Simon, to appear in "Journal of Computational and
 Graphical Statistics". The idea is as follows: it is possible to quickly
 check whether the solution to the graphical lasso problem will be block
 diagonal, for a given value of the tuning parameter. If so, then one can
 simply apply the graphical lasso algorithm to each block separately,
 leading to massive speed improvements.
 }
 
 \value{
 A list with components
 \item{w}{Estimated covariance matrix}
 \item{wi}{Estimated inverse covariance matrix}
 \item{loglik}{Value of maximized log-likelihodo+penalty}
 \item{errflag}{Memory allocation error flag: 0 means no error; !=0 means
  memory allocation error - no output returned}
 \item{approx}{Value of input argument approx}
 \item{del}{Change in parameter value at convergence}
 \item{niter}{Number of iterations of outer loop used by algorithm}
 }

\references{
 Bollhöfer, M., Eftekhari, A., Scheidegger, S. and Schenk, O., 2019. Large-scale sparse inverse covariance   matrix estimation. SIAM Journal on Scientific Computing, 41(1), pp.A380-A401.\url{https://epubs.siam.org/doi/abs/10.   1137/17M1147615?journalCode=sjoce3}

Eftekhari, A., Bollhöfer, M. and Schenk, O., 2018, November. Distributed memory sparse inverse covariance      matrix estimation on high-performance computing architectures. In SC18: International Conference for High   Performance Computing, Networking, Storage and Analysis (pp. 253-264). IEEE.  \url{https://dl.acm.org/doi/10.5555/       3291656.3291683}

Eftekhari A, Pasadakis D, Bollhöfer M, Scheidegger S, Schenk O (2021). Block-Enhanced PrecisionMatrix          Estimation for Large-Scale Datasets. Journal of Computational Science, p. 101389. \url{https://www.sciencedirect.com/   science/article/pii/S1877750321000776}
}

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

 out<-SQUIC(data,lambda)
 }


 \keyword{covariance matrices}
 \keyword{multivariate}
 \keyword{models}
 \keyword{graphs}

# SQUIC_R

### Sparse Quadratic Inverse Covariance Estimation

This is the SQUIC algorithm, made available as a Python package. SQUIC tackles the statistical problem of estimating large sparse inverse covariance matrices. This estimation poses an ubiquitous problem that arises in many applications e.g. coming from the fields mathematical finance, geology and health. SQUIC belongs to the class of second-order L1-regularized Gaussian maximum likelihood methods and is especially suitable for high-dimensional datasets with limited number of samples. For further details please see the listed references.

### Installation

1) Download the SQUIC shared library from www.gitlab.ci.inf.usi.ch/SQUIC/libSQUIC, and follow its README instructions.

2) Run the following command to install the library:

```angular2
library(devtools)  
install_github("www.gitlab.ci.inf.usi.ch/SQUIC/SQUIC_R")
```

### Example

To run a simple example : 

```angular2

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
```

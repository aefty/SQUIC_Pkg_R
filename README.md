# SQUIC R Interface Package

SQUIC is a second-order, L1-regularized maximum likelihood method for performant large-scale sparse precision matrix estimation. This repository contains the source code for the R interface of SQUIC. 

## Installation

Step 1: Download the shared library libSQUIC from www.gitlab.ci.inf.usi.ch/SQUIC/libSQUIC, and follow its README instructions. The default and recommended location for libSQUIC is the home directory, i.e., ``~/``.

Step 2: Run the following command to install the library:
```angular2
library(devtools) 
install_github("www.gitlab.ci.inf.usi.ch/SQUIC/SQUIC_R")
```
_Note: The ``devtools`` package can be install via the command ``install.packages("devtools")``._

Step 3: Load the SQUIC package:
```angular2
library(SQUIC)  
```
For further details type ``help(SQUIC)`` in the R command line.

_Note: The number of threads used by SQUIC can be defined by setting the enviroment variable OMP_NUM_THREADS (e.g., ``base> export OMP_NUM_THREADS=12``). This may require a restart of the session)._

## Example

To run a simple example : 

```angular2
library(SQUIC)

p = 1024
n = 100
lambda = .4

# generate a tridiagonal matrix
iC_star = Matrix::bandSparse(p, p, (-1):1, list(rep(-.5, p-1), rep(1.25,   p), rep(-.5, p-1)));

# generate the data
z    = replicate(n,rnorm(p));
iC_L = chol(iC_star);
data = matrix(solve(iC_L,z),p,n);

# Run SQUIC
out <- SQUIC(data,lambda)
```

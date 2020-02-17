# conv-rel-comb

Matlab code for reproducing the experimental results of the paper [Combinatorial Penalties:
Which structures are preserved by convex relaxations?](https://arxiv.org/abs/1710.06273), by Marwa El Halabi, Francis Bach, and Volkan Cevher.

## Usage:

To reproduce the figures in the paper run the script `AdaptiveTest_plot.m` (with `corr_level = 0` for uncorrelated setting, and `corr_level = 0.5` for correlated setting).

To reproduce the results call the function:

`AdaptiveTest(corr_level,n,tol,nits)`

with the following parameters (see paper for details):
* `corr_level` correlation level of measurements
*  `n` number of measurements
* `tol` threshold value - estimated support is the set of coefficients above this value (set to `1e-8`)
* `nruns` number of runs

**Dependencies**: CVX [GB14], and SPAMS toolbox [M17].

## References:

* [GB14] M. Grant and S. Boyd. CVX: Matlab software for disciplined convex programming, version 2.1. http://cvxr.com/cvx/, Mar. 2014.
* [M17]  Julien Mairal. SPAMS: Sparse Modeling Software, version 2.6. http://spams-devel.gforge.inria.fr/, retrieved Feb. 2017.



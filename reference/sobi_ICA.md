# SOBI ICA

Internal function for running SOBI ICA on an `eeg_epochs` object

## Usage

``` r
sobi_ICA(data, maxiter, tol, pca, centre, verbose = TRUE, whitening = "pca")
```

## Arguments

- data:

  Data to be ICAed.

- maxiter:

  Maximum number of iterations of the joint diagonalization

- tol:

  convergence tolerance.

- pca:

  Number of PCA components.

- centre:

  Mean centre signals.

- verbose:

  Print informative messages.

- whitening:

  Defaults to pca, options are pca or zca

## Author

A. Belouchrani and A. Cichocki. Adapted to R by Matt Craddock
<matt@mattcraddock.com>

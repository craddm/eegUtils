# Calculate coefficients and standard errors

Calculate coefficients and standard errors

## Usage

``` r
get_lm(mdf, x, inverted_dm, resid_df, robust = FALSE)
```

## Arguments

- mdf:

  model matrix

- x:

  data to be modelled (matrix of trials x channels)

- inverted_dm:

  inverted matrix - unscaled covariance

- resid_df:

  residual degrees of freedom for the model

- robust:

  Use heteroskedasticity-consistent covariance (HC3).

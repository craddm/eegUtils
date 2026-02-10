# Fit a linear model to EEG data

Fits a linear model to each timepoint separately for each electrode.

## Usage

``` r
fit_glm(formula, data, ...)

# S3 method for class 'eeg_epochs'
fit_glm(formula, data, baseline = NULL, ...)
```

## Arguments

- formula:

  An object of class `formula`. Right-hand side A regression formula for
  a GLM. See ?formula and notes on use below

- data:

  An `eegUtils` object.

- ...:

  Any other arguments passed to (LM/GLM)

- baseline:

  Numeric vector of length 2 specifying time period to be used as a
  baseline.

## Methods (by class)

- `fit_glm(eeg_epochs)`: GLM fitting for `eeg_epochs`

## Notes On Usage

The `fit_glm` function will fit a linear model to each timepoint for
each electrode in the input dataset.

The formula is a standard R formula. Specify only the right-hand side of
the formula i.e. any predictors.

The function allows flexible fitting of a baseline covariate,
recognising the term `baseline` in the specified formula.

## Author

Matt Craddock, <matt@mattcraddock.com>

## Examples

``` r
fit_glm(~epoch_labels, demo_spatial)
#> EEG Linear Model
#> 
#> Formula: ~ epoch_labels 
#> 
#> Number of fitted channels    :    23 
#> Channel names            : Fp1 F3 F7 C3 T7 P3 P7 Oz Pz Fpz Fp2 Fz F4 F8 Cz C4 T8 P4 P8 EXG1 EXG2 EXG3 EXG4 
#> Epoch limits         : -0.301 - 0.488 seconds
fit_glm(~epoch_labels + baseline, demo_spatial, baseline = c(-.1, 0))
#> EEG Linear Model
#> 
#> Formula: ~ epoch_labels + baseline 
#> 
#> Number of fitted channels    :    23 
#> Channel names            : Fp1 F3 F7 C3 T7 P3 P7 Oz Pz Fpz Fp2 Fz F4 F8 Cz C4 T8 P4 P8 EXG1 EXG2 EXG3 EXG4 
#> Epoch limits         : -0.301 - 0.488 seconds
```

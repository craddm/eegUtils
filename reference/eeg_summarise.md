# Calculate simple summary statistics for `eeg_*` objects

Calculate the timepoint-by-timepoint mean, standard deviation, standard
error, or variance `eeg_epochs` objects.

## Usage

``` r
eeg_summarise(data, ...)

# S3 method for class 'eeg_epochs'
eeg_summarise(
  data,
  statistic = c("sem", "mean", "sd", "var"),
  conditions = NULL,
  time_lim = NULL,
  ...
)
```

## Arguments

- data:

  An `eegUtils` object.

- ...:

  Various arguments passed to specific functions

- statistic:

  The statistic to calculate at each timepoint. Defaults to "sem"

- conditions:

  Conditions to group the data by.

- time_lim:

  Timepoint(s) to summarise. Can be a range, for which a summary
  statistic will be provided for each timepoint, or a list of individual
  times. If none is supplied, the function will calculate a summary for
  every timepoint.

## Value

A tibble

## Methods (by class)

- `eeg_summarise(eeg_epochs)`: Calculate summary statistics for
  `eeg_epochs` objects

## Examples

``` r
eeg_summarise(demo_spatial, statistic = "sem")
#> # A tibble: 102 × 24
#>      time sem_Fp1 sem_Fp2 sem_Fz sem_Cz sem_Pz sem_P3 sem_P4 sem_C3 sem_C4
#>     <dbl>   <dbl>   <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
#>  1 -0.301  0.0861  0.0769 0.0800 0.0933  0.108 0.0782 0.0879 0.0609 0.0775
#>  2 -0.293  0.0830  0.0763 0.0754 0.0981  0.105 0.0749 0.0837 0.0591 0.0733
#>  3 -0.285  0.0862  0.0733 0.0751 0.100   0.101 0.0737 0.0782 0.0598 0.0688
#>  4 -0.277  0.0907  0.0798 0.0755 0.0989  0.100 0.0730 0.0742 0.0664 0.0675
#>  5 -0.270  0.0885  0.0769 0.0743 0.0913  0.107 0.0726 0.0701 0.0647 0.0660
#>  6 -0.262  0.0857  0.0824 0.0713 0.0894  0.111 0.0795 0.0693 0.0610 0.0679
#>  7 -0.254  0.0860  0.0813 0.0763 0.0916  0.106 0.0859 0.0732 0.0598 0.0760
#>  8 -0.246  0.0877  0.0733 0.0823 0.0914  0.103 0.0823 0.0760 0.0563 0.0833
#>  9 -0.238  0.0835  0.0713 0.0793 0.0877  0.104 0.0730 0.0762 0.0563 0.0867
#> 10 -0.230  0.0847  0.0773 0.0704 0.0821  0.107 0.0739 0.0771 0.0606 0.0818
#> # ℹ 92 more rows
#> # ℹ 14 more variables: sem_P7 <dbl>, sem_P8 <dbl>, sem_F3 <dbl>, sem_F4 <dbl>,
#> #   sem_T7 <dbl>, sem_T8 <dbl>, sem_F7 <dbl>, sem_F8 <dbl>, sem_Oz <dbl>,
#> #   sem_Fpz <dbl>, sem_EXG1 <dbl>, sem_EXG2 <dbl>, sem_EXG3 <dbl>,
#> #   sem_EXG4 <dbl>
```

# Function to create an S3 object of class "eeg_stats".

Function to create an S3 object of class "eeg_stats".

## Usage

``` r
eeg_stats(statistic, chan_info, pvals, timings, method)
```

## Arguments

- statistic:

  Calculated statistic (e.g. t-statistic)

- chan_info:

  String of character names for electrodes.

- pvals:

  calculated p-values for that statistic

- timings:

  Unique timepoints remaining in the data.

- method:

  Type of statistical test

## Author

Matt Craddock <matt@mattcraddock.com>

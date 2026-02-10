# FASTER EEG artefact rejection

An implementation of the FASTER artefact rejection method for EEG by
Nolan, Whelan & Reilly (2010) FASTER: Fully Automated Statistical
Thresholding for EEG artifact Rejection. J Neurosci Methods.

## Usage

``` r
ar_FASTER(data, ...)

# S3 method for class 'eeg_epochs'
ar_FASTER(
  data,
  exclude = NULL,
  test_chans = TRUE,
  test_epochs = TRUE,
  test_cine = TRUE,
  ...
)

# S3 method for class 'eeg_group'
ar_FASTER(
  data,
  exclude = NULL,
  test_chans = TRUE,
  test_epochs = TRUE,
  test_cine = TRUE,
  EOG = NULL,
  ...
)

eeg_FASTER(data, ...)
```

## Arguments

- data:

  An object of class `eeg_epochs`

- ...:

  Parameters passed to FASTER

- exclude:

  Channels to be ignored by FASTER.

- test_chans:

  Logical. Run tests of global channel statistics

- test_epochs:

  Logical. Run tests of globally bad epochs.

- test_cine:

  Logical. Run tests for locally bad channels within epochs.

- EOG:

  names of EOG channels to be used when computed maximum EOG values.

## Value

An `eeg_epochs` object with artefact correction applied.

## Methods (by class)

- `ar_FASTER(eeg_epochs)`: Run FASTER on `eeg_epochs`

- `ar_FASTER(eeg_group)`: Run FASTER on `eeg_group` objects

## References

Nolan, Whelan & Reilly (2010). FASTER: Fully Automated Statistical
Thresholding for EEG artifact Rejection. J Neurosci Methods.

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
ar_FASTER(demo_epochs)
#> Globally bad channels: 
#> Globally bad epochs: 1 26 65
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 77 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
```

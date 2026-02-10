# Chan_info checker

Performs several checks on the structure of channel info: 1) Checks that
"electrode" is character, not factor. 2) rounds any numeric values to 2
decimal places. 3) Checks for any missing channels in the chan_info if
signal names are supplied; populates them with NA if it finds any.

## Usage

``` r
validate_channels(chan_info, sig_names = NULL)
```

## Arguments

- chan_info:

  A channel info structure

- sig_names:

  signal names from eegUtils signals

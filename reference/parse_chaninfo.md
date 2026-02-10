# Parse channel info from an EEGLAB set file

Internal function to convert EEGLAB chan_info to eegUtils style

## Usage

``` r
parse_chaninfo(chan_info, drop = FALSE)
```

## Arguments

- chan_info:

  Channel info list from an EEGLAB set file

- drop:

  If there are additional columns, remove all columns except electrode
  if TRUE, or just unexpected columns if FALSE.

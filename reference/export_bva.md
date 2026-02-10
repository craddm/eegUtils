# Export continuous data in Brain Vision Analyzer format

Export continuous EEG data in Brain Vision Analyzer format. This is one
of the recommended formats for BIDS
<https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html>

## Usage

``` r
export_bva(.data, filename, orientation, verbose = TRUE)

# S3 method for class 'eeg_data'
export_bva(.data, filename, orientation = "VECTORIZED", verbose = TRUE)
```

## Arguments

- .data:

  `eeg_data` object to be exported.

- filename:

  String giving filename to export to. File extensions will be removed
  when supplied.

- orientation:

  VECTORIZED or MULTIPLEXED. This relates to the way the data is stored
  in the binary file. VECTORIZED is the default and recommended.

- verbose:

  print informative messages to console

## Methods (by class)

- `export_bva(eeg_data)`: Method for `eeg_data`

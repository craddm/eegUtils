# Internal function for running SSD algorithm

Internal function for running SSD algorithm

## Usage

``` r
run_SSD(data, sig_range, noise_range, RESS = FALSE, verbose = TRUE, order = 2)
```

## Arguments

- data:

  `eeg_epochs` object to be decomposed

- sig_range:

  Frequency range of the signal of interest

- noise_range:

  Frequency range of the noise

- RESS:

  Run RESS rather than SSD. Defaults to FALSE.

- verbose:

  Informative messages in consoles. Defaults to TRUE.

- order:

  filter order for IIR filters

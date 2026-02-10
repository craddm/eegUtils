# Calculate inter-trial coherence

Calculates inter-trial coherence (ITC), a measure of phase consistency
across single trial data. Input data must be provided as complex Fourier
coefficients within an `eeg_tfr` object

## Usage

``` r
compute_itc(data)
```

## Arguments

- data:

  An `eeg_tfr` object

## Value

An `eeg_tfr` object

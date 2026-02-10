# EEG decomposition viewer

A Shiny viewer for Independent Component Analysis or Spatio-spectral
Decomposition/RESS components that provides an interface for looking at
topographies, timecourses, and power spectral densities of all or
individual components. Can be used to select and reject artefactual
components.

## Usage

``` r
view_ica(data)
```

## Arguments

- data:

  An `eeg_ICA` object

## Value

A list consisting (optionally) of

- A character vector of components marked for rejection

- A character vector of components marked to be kept

- An `eeg_epochs` object reconstructed from the `eeg_ICA` object, with
  components marked for rejection removed.

## Author

Matt Craddock <matt@mattcraddock.com>

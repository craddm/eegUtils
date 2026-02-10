# Channel interpolation

Interpolate EEG channels using a spherical spline (Perrin et al., 1989;
1990). The data must have channel locations attached.

## Usage

``` r
interp_elecs(data, bad_elecs, ...)

# S3 method for class 'eeg_data'
interp_elecs(data, bad_elecs, ...)
```

## Arguments

- data:

  Data as a `eeg_data` or `eeg_epochs` object.

- bad_elecs:

  Name(s) of electrode(s) to interpolate.

- ...:

  Other parameters passed to the functions.

## Methods (by class)

- `interp_elecs(eeg_data)`: Interpolate EEG channel(s)

## References

- Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F. (1989).
  Spherical splines for scalp potential and current density mapping.
  Electroencephalography and Clinical Neurophysiology, 72, 184-187

- Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F. (1990).
  Corrigenda EEG 02274. Electroencephalography and Clinical
  Neurophysiology, 76, 565

## Author

Matt Craddock <matt@mattcraddock.com>

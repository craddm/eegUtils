# Normalise Morlet wavelets

Normalise the FFT'd Morlet wavelet family as recommended by Mike X
Cohen, dividing each wavelet by its absolute maximum. This should result
in each frequency being passed at unit amplitude, and the resulting
convolution with the signal should return units approximately on the
original scale (i.e. uV^2 / Hz). an alternative would be Frobenius norm?

## Usage

``` r
wavelet_norm(mf_zp, n_freq)
```

## Arguments

- mf_zp:

  A zero-padded, FFT'd morlet wavelet family

- n_freq:

  Number of frequencies

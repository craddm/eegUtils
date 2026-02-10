# A guide to eegUtils data structures

## eegUtils

`eegUtils` uses S3 objects to store EEG data and associated information
such as channel locations. Using different object classes for data
structured in different ways ensures that the various plotting functions
work consistently across different types of EEG data. For example, there
are different classes for epoched (`eeg_epochs`) and continuous data
(`eeg_data`), and for time-frequency representations of data
(`eeg_tfr`).

## `eeg_data` objects

`eeg_data` objects are the base class used for continuous data. When raw
data is imported, the output is this class. This class is a list
constituting the following entries:

- `signals`
  - A data frame containing the actual EEG data in wide format. Each
    column is data from a single electrode, each row from a single
    timepoint.
- `srate`
  - A single numeric value giving the sampling rate of the data in Hz
- `events`
  - A data.frame/tibble with 3 columns describing the events recorded in
    the data.
    - event_onset (in samples) relative to *recording onset*
    - event_time (in seconds) relative to *recording onset*
    - event_type (typically integer, others possible)
- `chan_info`
  - Often NULL on import, since BDF and CNT do not contain (usable)
    electrode locations.
  - A data.frame/tibble containing channel location information can be
    added. Currently similar to EEGLAB style, with the following columns
    (may change):
    - `electrode` - electrode names
    - `radius` - Spherical co-ordinates (Radius is typically normalized
      to 1)
    - `theta` - Spherical co-ordinates (theta)
    - `phi` - Spherical co-ordinates (theta)
    - `cart_x` - Cartesian 3D coordinates
    - `cart_y` - Cartesian 3D coordinates
    - `cart_z` - Cartesian 3D coordinates
    - `x` - 2D Stereographic projection of the spherical coordinates
    - `y` - 2D Stereographic projection of the spherical coordinates

    ``` r
    head(eegUtils:::electrodeLocs)
    #> # A tibble: 6 × 9
    #>   electrode radius theta   phi cart_x cart_y cart_z      x     y
    #>   <fct>      <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>
    #> 1 LPA            1  -119    13 -86.1   -20.0 -48.0  -116.  -26.8
    #> 2 RPA            1   119   -13  85.8   -20.0 -48.0   116.  -26.8
    #> 3 Nz             1   115    90   0.01   86.8 -40.0     0   115  
    #> 4 Fp1            1   -94   -71 -29.4    83.9  -6.99  -30.6  88.9
    #> 5 Fpz            1    91    90   0.11   88.2  -1.71    0    91  
    #> 6 Fp2            1    94    71  29.9    84.9  -7.08   30.6  88.9
    ```
- `timings`
  - A data.frame containing a description of each row in time (s) and
    sample numbers (samples)
- `epochs`
  - A data.frame containing information about the epochs in the data,
    more relevant for epoched data.
- `reference`
  - NA on loading, in the absence of a recorded reference channel.
  - Once data is referenced, is a list consisting of two entries.
    - `ref_chans` - Labels for channels used to calculate the reference
      data. Can also be “average”.
    - `excluded` - Labels for any channels excluded from the reference
      data.

## eeg_epochs

`eeg_epochs` objects share the same overall structure with `eeg_data`
objects, but some of the internals currently differ, as described below.

- `events`
  - The events table has two additional columns, `epoch` and `time`.
    - `epoch` gives the epoch number to which a given event belongs
    - `time` gives the time point at which the event occurs relative to
      the *epoch onset*
    - `event_time` still gives the time point at which the event occurs
      relative to the *recording onset*

  ``` r
  events(demo_epochs)
  #> # A tibble: 80 × 5
  #>    event_onset event_time event_type epoch  time
  #>          <dbl>      <dbl>      <dbl> <dbl> <dbl>
  #>  1        4128       8.06        208     1     0
  #>  2        7037      13.7         213     2     0
  #>  3       10043      19.6         215     3     0
  #>  4       12928      25.2         213     4     0
  #>  5       15868      31.0         207     5     0
  #>  6       18777      36.7         207     6     0
  #>  7       21578      42.1         213     7     0
  #>  8       24554      48.0         213     8     0
  #>  9       27379      53.5         222     9     0
  #> 10       30306      59.2         208    10     0
  #> # ℹ 70 more rows
  ```
- `timings`
  - The timings table has one additional column, `epoch`.
    - `epoch` gives the epoch number to which a given datapoint belongs
    - `sample` still uniquely identifies each datapoint
    - `time` now gives the time relative to the zero-point of the epoch,
      i.e. the *event* on which the epoch is centred.

  ``` r
  head(demo_epochs$timings)
  #> # A tibble: 6 × 3
  #>   epoch sample   time
  #>   <dbl>  <dbl>  <dbl>
  #> 1     1   4028 -0.197
  #> 2     1   4032 -0.189
  #> 3     1   4036 -0.182
  #> 4     1   4040 -0.174
  #> 5     1   4044 -0.166
  #> 6     1   4048 -0.158
  ```
- `epochs`
  - A data.frame containing information about the participant, and
    labels applied to any epochs

  ``` r
  epochs(demo_epochs)
  #> # A tibble: 80 × 4
  #>    epoch recording epoch_label participant_id
  #>    <dbl> <lgl>     <lgl>       <chr>         
  #>  1     1 NA        NA          001           
  #>  2     2 NA        NA          001           
  #>  3     3 NA        NA          001           
  #>  4     4 NA        NA          001           
  #>  5     5 NA        NA          001           
  #>  6     6 NA        NA          001           
  #>  7     7 NA        NA          001           
  #>  8     8 NA        NA          001           
  #>  9     9 NA        NA          001           
  #> 10    10 NA        NA          001           
  #> # ℹ 70 more rows
  ```

## eeg_tfr

`eeg_tfr` objects hold time-frequency representations of `eeg_epochs`
objects.

- `signals`
  - a 4D matrix - epoch by time by electrode by frequency
- `dimensions`
  - keeps track of which matrix dimension corresponds to which property.

## eeg_ICA

`eeg_ICA` objects contain the results of either an ICA or an SSD
decomposition applied to an `eeg_epochs` object.

- `mixing_matrix`
  - The weights that are used to convert the source data into the
    original data
- `unmixing_matrix`
  - The weights that are used to convert the original data into source
    data
- `signals` - individual component activations

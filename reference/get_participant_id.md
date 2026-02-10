# Query and set elements of the `epochs` metadata structures

These functions are used to query and set the `participant_id` and
`recording` fields of an `epochs` structure within any `eegUtils`
object. Note that the two `set_*` functions do not operate on
`eeg_group` objects.

## Usage

``` r
get_participant_id(data)

get_recording(data)

set_participant_id(data, participant_id)

set_recording(data, recording)
```

## Arguments

- data:

  An `eegUtils` object.

- participant_id:

  A character vector giving the name to use for the `participant_id` for
  a particular object.

- recording:

  A character vector giving the name to use for the EEG recording.

## Value

- `get_participant_id()` returns a character vector containing the
  `participant_id`

- `get_recording()` returns a character vector containing the
  `recording` from values of `recording` from the `eegUtils` object

- `set_participant_id()` returns the full `eegUtils` object with the
  `participant_id` field modified in the `epochs` structure

- `set_recording()` returns the full `eegUtils`object with the
  `recording` field modified in the `epochs` structure

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
get_participant_id(demo_epochs)
#> [1] "001"
get_recording(demo_epochs)
#> [1] NA
set_participant_id(demo_epochs, "002")
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
set_recording(demo_epochs, "test_recording")
#> Epoched EEG data
#> 
#> Number of channels   : 11 
#> Number of epochs : 80 
#> Epoch limits     : -0.197 - 0.451 seconds
#> Electrode names      : A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
#> Sampling rate        : 128  Hz
#> Reference        : average 
```

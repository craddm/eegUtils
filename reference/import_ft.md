# Import Fieldtrip files

Fieldtrip is a Matlab package for EEG/MEG processing and analysis.

## Usage

``` r
import_ft(file_name, participant_id = NULL, recording = NULL, verbose = TRUE)
```

## Arguments

- file_name:

  Name of file to be imported.

- participant_id:

  Identifier for the participant.

- recording:

  Name of the recording. By default, the filename will be used.

- verbose:

  Informative messages printed to console. Defaults to TRUE.

## Value

An object of class `eeg_data`, `eeg_epochs`, or `eeg_tfr`, depending on
the type of input data.

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
if (FALSE) import_ft("fieldtrip_test.mat") # \dontrun{}
```

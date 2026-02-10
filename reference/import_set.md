# Load `EEGLAB` .set files

Load `EEGLAB` .set files and convert them to `eeg_epochs` objects.
Supports import of files saved in either Matlab v6.5 or Matlab v7.3
formats. Currently, any ICA weights or decompositions are discarded.

## Usage

``` r
import_set(
  file_name,
  df_out = FALSE,
  participant_id = NULL,
  recording = NULL,
  drop_custom = FALSE,
  verbose = TRUE
)
```

## Arguments

- file_name:

  Filename (and path if not in present working directory)

- df_out:

  Defaults to FALSE - outputs an object of class `eeg_epochs`. Set to
  TRUE for a normal data frame.

- participant_id:

  Character vector. By default, the filename will be used as the id of
  the participant.

- recording:

  Character vector. By default, the filename will be used as the name of
  the recording.

- drop_custom:

  Drop custom event fields. TRUE by default.

- verbose:

  Print informative messages. TRUE by default.

## Value

An object of class `eeg_epochs`

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
if (FALSE) import_set("your_data.set") # \dontrun{}
```

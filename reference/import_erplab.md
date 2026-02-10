# Import from ERPLAB .erp files

ERPLAB is a toolbox for Event Related Potential analysis associated with
EEGLAB. It exports files in Matlab v6.5 format, either with an .erp or
.mat extension. Note that ERPLAB files can contain data from multiple
subjects, individual subjects, or averaged over multiple subjects. The
files do not distinguish between these possibilities, so users will need
to apply their own knowledge of the contents to import them correctly.

## Usage

``` r
import_erplab(
  file_name,
  df_out = FALSE,
  participant_id = NULL,
  recording = NULL
)
```

## Arguments

- file_name:

  Filename (and path if not in present working directory)

- df_out:

  Defaults to FALSE - outputs an object of class `eeg_epochs`. Set to
  TRUE to output a data frame.

- participant_id:

  By default, the filename will be used as the id of the participant.

- recording:

  By default, the filename will be used as the name of the recording.

## Value

An object of class `eeg_epochs` or a data.frame.

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
if (FALSE) import_set("your_data.set") # \dontrun{}
```

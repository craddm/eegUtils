# Function for reading raw data.

Currently BDF/EDF, 32-bit .CNT, and Brain Vision Analyzer files are
supported. Filetype is determined by the file extension.The `edfReader`
package is used to load BDF/EDF files, whereas custom code is used for
.CNT and BVA files. The function creates an `eeg_data` structure for
subsequent use.

## Usage

``` r
import_raw(
  file_name,
  file_path = NULL,
  recording = NULL,
  participant_id = NA,
  fast_bdf = TRUE
)
```

## Arguments

- file_name:

  File to import. Should include file extension.

- file_path:

  Path to file name, if not included in filename.

- recording:

  Name of the recording. By default, the filename will be used.

- participant_id:

  Identifier for the participant. Defaults to NA.

- fast_bdf:

  New, faster method for loading BDF files. Experimental.

## Value

An object of class `eeg_data`

## Author

Matt Craddock, <matt@mattcraddock.com>

## Examples

``` r
if (FALSE) { # \dontrun{
import_raw("test_bdf.bdf")
} # }
```

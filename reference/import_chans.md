# Import channel locations from various file formats

Currently only ASA `.elc` format with Cartesian x-y-z coordinates is
supported.

## Usage

``` r
import_chans(file_name, format = "spherical", file_format = "auto")
```

## Arguments

- file_name:

  Name and full path of file to be loaded.

- format:

  If the file is not `.elc` format, "spherical", "geographic". Default
  is "spherical".

- file_format:

  Default is `auto`, which will use the file extension to determine file
  format. Other options include `ced`, `besa`, `elp`, `elc`

## Value

A `tibble` containing electrode names and locations in several different
coordinate systems.

## Author

Matt Craddock <matt@mattcraddock.com>

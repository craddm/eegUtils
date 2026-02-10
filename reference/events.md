# Modify events structure

Get or set the values in the `events` structure of an eegUtils object.

## Usage

``` r
events(.data)

events(.data) <- value

# S3 method for class 'eeg_epochs'
events(.data) <- value

# S3 method for class 'eeg_data'
events(.data) <- value
```

## Arguments

- .data:

  `eegUtils` object to view

- value:

  Value to replace `events` structure with.

## See also

Other event handlers:
[`list_epochs()`](https://craddm.github.io/eegUtils/reference/list_epochs.md),
[`list_events()`](https://craddm.github.io/eegUtils/reference/list_events.md),
[`tag_events()`](https://craddm.github.io/eegUtils/reference/tag_events.md)

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

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
events(demo_epochs) <- mutate(events(demo_epochs),
 sf = dplyr::case_when(
         event_type %% 2 == 0 ~ "HSF",
         event_type %% 2 == 1 ~ "LSF",
 ))
events(demo_epochs)
#> # A tibble: 80 × 6
#>    event_onset event_time event_type epoch  time sf   
#>          <dbl>      <dbl>      <dbl> <dbl> <dbl> <chr>
#>  1        4128       8.06        208     1     0 HSF  
#>  2        7037      13.7         213     2     0 LSF  
#>  3       10043      19.6         215     3     0 LSF  
#>  4       12928      25.2         213     4     0 LSF  
#>  5       15868      31.0         207     5     0 LSF  
#>  6       18777      36.7         207     6     0 LSF  
#>  7       21578      42.1         213     7     0 LSF  
#>  8       24554      48.0         213     8     0 LSF  
#>  9       27379      53.5         222     9     0 HSF  
#> 10       30306      59.2         208    10     0 HSF  
#> # ℹ 70 more rows
```

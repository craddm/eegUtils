# Modify the epochs structure

Get or set the epochs structure of an `eegUtils` object.

## Usage

``` r
epochs(data)

epochs(data) <- value
```

## Arguments

- data:

  `eegUtils` object to view

- value:

  Structure to replace `epochs` structure with.

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
  epochs(demo_spatial)
#> # A tibble: 80 × 5
#>    epoch participant_id recording         event_type epoch_labels
#>    <dbl> <chr>          <chr>                  <dbl> <chr>       
#>  1     1 ""             Matt-task-spatcue        122 valid_right 
#>  2     2 ""             Matt-task-spatcue        120 valid_left  
#>  3     3 ""             Matt-task-spatcue        122 valid_right 
#>  4     4 ""             Matt-task-spatcue        122 valid_right 
#>  5     5 ""             Matt-task-spatcue        120 valid_left  
#>  6     6 ""             Matt-task-spatcue        122 valid_right 
#>  7     7 ""             Matt-task-spatcue        120 valid_left  
#>  8     8 ""             Matt-task-spatcue        120 valid_left  
#>  9     9 ""             Matt-task-spatcue        122 valid_right 
#> 10    10 ""             Matt-task-spatcue        122 valid_right 
#> # ℹ 70 more rows
```

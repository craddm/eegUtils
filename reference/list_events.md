# List events

List trigger types and any labels found in an `eeg_data` object.

## Usage

``` r
list_events(data)
```

## Arguments

- data:

  An object of class `eeg_data`

## See also

[`tag_events()`](https://craddm.github.io/eegUtils/reference/tag_events.md)
and
[`list_epochs()`](https://craddm.github.io/eegUtils/reference/list_epochs.md)

Other event handlers:
[`events()`](https://craddm.github.io/eegUtils/reference/events.md),
[`list_epochs()`](https://craddm.github.io/eegUtils/reference/list_epochs.md),
[`tag_events()`](https://craddm.github.io/eegUtils/reference/tag_events.md)

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
list_events(demo_epochs)
#>   event_type
#> 1        208
#> 2        213
#> 3        215
#> 4        207
#> 5        222
#> 6        219
```

# Plot Global Field Power of EEG Signals

Global Field Power (Lehmann & Skrandies, 1980) is a way to quantify the
amount of activity in the overall electroencephalographic signal at a
given timepoint. It corresponds to the spatial standard deviation.

## Usage

``` r
plot_gfp(data, cols = NULL, keep_trials = FALSE)
```

## Arguments

- data:

  An `eeg_epochs` object

- cols:

  Condition columns from the `epochs` metadata to calculate GFP
  separately for different conditions.

- keep_trials:

  Calculate GFP for each epoch separately, then average over epochs.

## Author

Matt Craddock <matt@mattcraddock.com>

## Examples

``` r
plot_gfp(demo_spatial)
#> Creating epochs based on combinations of variables: participant_id 

plot_gfp(demo_spatial, keep_trials = TRUE)

plot_gfp(demo_spatial, cols = "epoch_labels")
#> Creating epochs based on combinations of variables: participant_id epoch_labels 

plot_gfp(demo_spatial, cols = "epoch_labels", keep_trials = TRUE)
```

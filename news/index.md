# Changelog

## eegUtils 0.8.0

#### Function changes

- [`plot_timecourse()`](https://craddm.github.io/eegUtils/reference/plot_timecourse.md)
  now requires facets and mappings to be explicitly stated during the
  call, rather than added afterwards. This allows it to use weighted
  averages when required.
- [`plot_timecourse.eeg_tfr()`](https://craddm.github.io/eegUtils/reference/plot_timecourse.md)
  now defaults to no baseline correction, and takes a `type` argument
  for baseline type to be specified.
- Added new
  [`plot_gfp()`](https://craddm.github.io/eegUtils/reference/plot_gfp.md)
  function for calculating and plotting Global Field Power.
- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  now supports plotting of multiple timepoints - pass a list of times to
  the `time_lim` argument.
- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  now supports custom fill titles through the `fill_title` argument, and
  automatically switches where necessary (e.g. now says “Power mV^2”).
  Fill titles are now also centred.
- [`browse_data.eeg_ICA()`](https://craddm.github.io/eegUtils/reference/browse_data.md)
  now provides the option to select components for rejection, and
  returns a character vector of selected components.
- Recoded
  [`browse_data()`](https://craddm.github.io/eegUtils/reference/browse_data.md)
  to use `bslib` for styling.
- [`view_ica()`](https://craddm.github.io/eegUtils/reference/view_ica.md)
  now allows you to select components for rejection, to double-click on
  topographies to inspect them individually, and to return cleaned data.
- [`eeg_average()`](https://craddm.github.io/eegUtils/reference/eeg_average.md)
  now supports averaging over conditions in `eeg_evoked` files
- [`eeg_average()`](https://craddm.github.io/eegUtils/reference/eeg_average.md)
  now records weights - the number of epochs that went into an average -
  and uses those in subsequent steps where possible for
  `eeg_epochs`/`eeg_evoked` / `eeg_tfr` objects.
- [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md)
  now has an argument `trim_edges` which allows users to switch off
  automatic removal of epoch edges after transformation. Defaults to
  TRUE.
- [`rm_baseline()`](https://craddm.github.io/eegUtils/reference/rm_baseline.md)
  now adds a record of the baseline period
- Changed logic of
  [`eeg_combine.eeg_epochs()`](https://craddm.github.io/eegUtils/reference/eeg_combine.md).
  Now checks for duplicate epochs directly using `epoch`,
  `participant_id` and `recording`, and only corrects when there are
  duplicates. This means `recording` and `participant_id` have to be
  matching across `eeg_epochs` objects for the timing correction to be
  applied.

#### Internal changes/bug fixes

- Recoded
  [`faster_epochs()`](https://craddm.github.io/eegUtils/reference/faster_epochs.md)
  to no longer use `data.table`.
- added `parse_cycles()` function for use during
  [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md)
- Using
  [`electrode_locations()`](https://craddm.github.io/eegUtils/reference/electrode_locations.md)
  on a `data.frame` would return a data frame with the electrode names
  in full upper case. Now returns with the electrodes in their original
  case.
- `eeg_evoked` objects should now contain reference information.
- `ar_for_ica` file started, moving some functions from
  `artefact_rejection.R`
- `select_times` prevented from converting single column data.frames to
  vectors after subsetting.
- `plot_timecourse.eeg_tfr` now correctly uses the `freq_range`
  parameter.
- Improvements to internal processing logic in
  [`view_ica()`](https://craddm.github.io/eegUtils/reference/view_ica.md)
  to improve performance.
- `import_raw(..., fast_bdf = TRUE)` will now discard Annotations rather
  than fail to import BDF files with Annotations.
- removed dependency on `Matrix` - now using `qr()$rank` in `run_ICA`
  rather than using
  [`Matrix::rankMatrix()`](https://rdrr.io/pkg/Matrix/man/rankMatrix.html)
  to determine rank of input signals.
- Added numerous `dropped_aes` variables to the custom `ggplot2` `stat_`
  functions for compatibility with `ggplot2 3.4.0`.
- Replaced `size` aesthetic with `linewidth` for compatibility with
  `ggplot2 3.4.0`.
- Lots of minor code style improvements
- Deprecated function `iir_filt` removed
- Refactored `run_ICA` and adopted faster SOBI methods
- Used internal legendre polynomial function and removed `pracma`
  dependency
- Used internal whitening methods in `sobi_ICA` and removed
  `whitenening` dependency

## eegUtils 0.7.0

#### Function changes

- Add `imax` method to
  [`run_ICA()`](https://craddm.github.io/eegUtils/reference/run_ICA.md).
  This allows use of the `infomax` ICA algorithm from the `infomax`
  package, which is a reimplementation of the extended-Infomax algorithm
  used in the `EEGLAB` Matlab toolbox.
- [`erp_scalp()`](https://craddm.github.io/eegUtils/reference/erp_scalp.md)
  and
  [`interactive_scalp()`](https://craddm.github.io/eegUtils/reference/interactive_scalp.md)
  should now appropriately use channel locations included in the data.
- More informative messages when using
  [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md).
- [`plot_tfr()`](https://craddm.github.io/eegUtils/reference/plot_tfr.md)
  now applies baseline correction on a single-trial basis where
  possible, which may show different results when using non-linear
  baseline correction (e.g. `divide` or `dB`)
- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  now allows you to provide multiple component numbers when plotting
  from an `eeg_ICA` object, and will automatically produce an
  appropriately facetted plot.
- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  now has a `k` parameter to control the smoothing when using
  `method = "gam"`.
- added additional `demo_spatial` data from a spatial cueing experiment.
- [`plot_difference()`](https://craddm.github.io/eegUtils/reference/plot_difference.md)
  function added for plotting ERP difference waves. Only currently
  handles two levels.
- Added `hanning` taper support for `compute_tfr`. Note that the scaling
  factors used for all `compute_tfr` calculations have been adjusted, so
  the exact numerical values returned will change. However, this is just
  a scaling factor - the relative distances between values remained
  unchanged.
- [`ar_FASTER()`](https://craddm.github.io/eegUtils/reference/ar_FASTER.md)
  has experimental support for `eeg_group` objects when those objects
  are `eeg_evoked` groups. It does not perform rejection but reports how
  many times each participants data breaks a threshold for a number of
  measures.
- [`import_raw()`](https://craddm.github.io/eegUtils/reference/import_raw.md)
  default `participant_id` is now changed to `NA` instead of
  `character(1)`, to promote better use with
  [`eeg_combine()`](https://craddm.github.io/eegUtils/reference/eeg_combine.md).
- [`eeg_combine()`](https://craddm.github.io/eegUtils/reference/eeg_combine.md)
  will now refuse to combine objects where `participant_id` is missing
  (i.e. is `NA`), and warn when combining objects with the previous
  default value ““. This is to prevent accidentally treating data from
  different participants as being from the same participant.
- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  now provides informative messages about the head radius used for
  plotting. The default calculation of `r` when using
  `interp_limit = "head"` has changed and should now be set at the
  outermost electrode’s position + a 10% buffer. Smaller `r` can be set
  manually.
- added `get_participant_id`, `set_participant_id`, `get_recording` and
  `set_recording` to interact with the `epochs` metadata in each
  `eegUtils` object.
- [`eeg_combine()`](https://craddm.github.io/eegUtils/reference/eeg_combine.md)
  handles
  [`eeg_tfr()`](https://craddm.github.io/eegUtils/reference/eeg_tfr.md)
  objects better, now returns an error when trying to combine
  single-trial data.
- [`plot_timecourse()`](https://craddm.github.io/eegUtils/reference/plot_timecourse.md)
  now handles
  [`eeg_tfr()`](https://craddm.github.io/eegUtils/reference/eeg_tfr.md)
  objects.

#### Internal changes / bug fixes

- [`erp_scalp()`](https://craddm.github.io/eegUtils/reference/erp_scalp.md)
  and
  [`interactive_scalp()`](https://craddm.github.io/eegUtils/reference/interactive_scalp.md)
  now use cleaner evaluation of the `colour` argument using `rlang`
- Fixed adding CSD as new reference when using
  [`compute_csd()`](https://craddm.github.io/eegUtils/reference/compute_csd.md)
- When using a scaling number of cycles in
  [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md),
  they will now also use the `spacing` parameter to determine `log` or
  `linear` scaling.
- Fixed bug with
  [`eeg_average()`](https://craddm.github.io/eegUtils/reference/eeg_average.md)
  used on `c("eeg_group", "eeg_tfr")` objects.
- Fixed bug with incorrect number of epochs calculated when epoching
  `eeg_data` objects if there were multiple “target” triggers appearing
  in an epoch.
- Fixed error with
  [`erp_raster()`](https://craddm.github.io/eegUtils/reference/erp_raster.md)
  when `anat_order == FALSE`
- new `epoch_queries` file for functions for setting and getting the
  `epochs` structure
- new `check_items` file for functions that check for consistency of
  various structures
- [`eeg_combine()`](https://craddm.github.io/eegUtils/reference/eeg_combine.md)
  with
  [`eeg_tfr()`](https://craddm.github.io/eegUtils/reference/eeg_tfr.md)
  no longer drops single dimensions, which was causing issues when there
  was only one channel or epoch in the object.
- `plot_timecourse.eeg_tfr` now correctly passes baseline period to
  [`rm_baseline()`](https://craddm.github.io/eegUtils/reference/rm_baseline.md)
- `filter.eeg_tfr()` was sometimes dropping single dimensions when using
  [`abind::asub()`](https://rdrr.io/pkg/abind/man/asub.html), now fixed

## eegUtils 0.6.3

#### Function changes

- Added log spaced frequencies to
  [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md),
  with the new `spacing` argument. `plot_tfr` automatically detects the
  spacing and plots the figure appropriately.
- Added `na.rm` option to
  [`erp_image()`](https://craddm.github.io/eegUtils/reference/erp_image.md)
  to either keep or plot trials with NA values due to smoothing. By
  default they’ll be removed.
- Added more informative messages for
  [`compute_psd()`](https://craddm.github.io/eegUtils/reference/compute_psd.md).

#### Internal changes / bug fixes

- Some minor documentation fixes.
- [`plot_tfr()`](https://craddm.github.io/eegUtils/reference/plot_tfr.md)
  error when selecting a specific frequency range fixed.
- Switched to new style of `vdiffr` tests
- fixed
  [`erp_image()`](https://craddm.github.io/eegUtils/reference/erp_image.md)
  smoothing over time instead of epochs.

## eegUtils 0.6.2

#### Function changes

- added support for `EEGLAB` .set files saved in newer Matlab file
  formats.
- changed first argument of
  [`eeg_filter()`](https://craddm.github.io/eegUtils/reference/eeg_filter.md)
  to `data` instead of `.data`
- added some more informative user messages for importing data and
  adding electrode locations.
- [`erp_image()`](https://craddm.github.io/eegUtils/reference/erp_image.md)
  now supports `eeg_tfr` objects.

#### Internal changes / bug fixes

- When combining three or more continuous `eeg_data` objects,
  [`eeg_combine()`](https://craddm.github.io/eegUtils/reference/eeg_combine.md)
  would substantially undercorrect the timing of events in the third
  file - this is now fixed.
- `groups` parameter for
  [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  now correctly passed for all types of object.
- `stat_tests.R` file removed, will be reimplemented elsewhere
- Long standing issues with import of channel locations from EEGLAB
  files hopefully fixed…
- [`rm_baseline()`](https://craddm.github.io/eegUtils/reference/rm_baseline.md)
  for `eeg_evoked` no longer uses `data.table`
- [`as.data.frame.eeg_evoked()`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_evoked.md)
  handles grouped data better.
- [`import_set()`](https://craddm.github.io/eegUtils/reference/import_set.md)
  now handles all EEGLAB formats better.

## eegUtils 0.6.1

#### Function changes

- `ar_eogcor` now has a `bipolarize` argument which can be set to false
  when the HEOG/VEOG channels are already bipolarized.
- added some new `ggplot2` based functions for topoplotting and adding
  contours
  - [`stat_scalpcontours()`](https://craddm.github.io/eegUtils/reference/stat_scalpcontours.md)
- [`geom_topo()`](https://craddm.github.io/eegUtils/reference/geom_topo.md)
  now has contours
- added errors with when attempting to use `compute_psd` or
  `compute_tfr` on `eeg_group` objects.
- [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md)
  now works better with `eeg_evoked` objects that contain multiple
  conditions.

#### Internal changes / bug fixes

- `eeg_reference` now handles multiple reference channels better on
  rereferencing
- `get_scalpmap` handles `eeg_ICA` components better when there are
  channels with no locations
- Travis-CI removed.
- `cart_to_spherical` coord flipping bug fixed (hopefully…)
- `filter` now converts to tibble internally and does not coerce
  `signals` to a vector when there is only one channel.
- added copyright info to `summary_contour` file
- `eeg_combine.eeg_evoked` made to behave more consistently when
  creating grouped data
- all [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  functions moved to `df_converters.r`

## eegUtils 0.6.0

IMPORTANT: There have been some changes to the logic of the
[`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
that may make their appearance quite different. Specifically, these
changes are to the way the underlying interpolation grid is calculated
and to how things like the diameter of the cartoon head is calculated.
These changes often lead to different minimum or maximum amplitudes
across the image, and thus changes in the appearance of the plot due to
different scales- don’t be alarmed!

#### Function changes

- `epoch_data` baseline correction now defaults to *no* baseline
  correction
- Added `filter` method for `eeg_tfr` objects
- [`fit_glm()`](https://craddm.github.io/eegUtils/reference/fit_glm.md)
  overhauled. Now far faster and allows specification of models using
  standard R formulae.
- New `eeg_lm` class introduced for output of
  [`fit_glm()`](https://craddm.github.io/eegUtils/reference/fit_glm.md).
- [`plot_butterfly.eeg_lm()`](https://craddm.github.io/eegUtils/reference/plot_butterfly.md)
  method added
- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  have been added for `eeg_lm` objects.
- [`view_ica()`](https://craddm.github.io/eegUtils/reference/view_ica.md)
  Shiny viewer for `eeg_ICA` and `eeg_decomp` objects added.
- [`view_artefacts()`](https://craddm.github.io/eegUtils/reference/view_artefacts.md)Shiny
  viewer for channel and epoch stats added.
- [`plot_timecourse()`](https://craddm.github.io/eegUtils/reference/plot_timecourse.md)
  now takes a mapping argument, which allows use of `ggplot2`
  [`aes()`](https://ggplot2.tidyverse.org/reference/aes.html) mappings
- [`eeg_average.eeg_tfr()`](https://craddm.github.io/eegUtils/reference/eeg_average.md)
  now follows behaviour of other
  [`eeg_average()`](https://craddm.github.io/eegUtils/reference/eeg_average.md)
  methods in respecting the `epochs` structure.  
- [`compute_itc()`](https://craddm.github.io/eegUtils/reference/compute_ITC.md)
  added for computing inter-trial coherence from `eeg_tfr` objects.
- `cols` added to `eeg_average.eeg_tfr`
- `eeg_combine.tfr_average()` added to handle pre-averaged `eeg_tfr`
  objects
- [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md)
  now allows non-constant number of cycles
- [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md)
  now uses a different scaling factor, so raw units should now be
  microvolts-squared.
- added
  [`import_erplab()`](https://craddm.github.io/eegUtils/reference/import_erplab.md)
  function
- [`plot_timecourse()`](https://craddm.github.io/eegUtils/reference/plot_timecourse.md)
  now allows CIs for `eeg_group` objects.

#### Internal changes / bug fixes

- [`plot_tfr()`](https://craddm.github.io/eegUtils/reference/plot_tfr.md)
  now always drops NA/NaN values and averages appropriately over
  electrodes and conditions.
- [`import_set()`](https://craddm.github.io/eegUtils/reference/import_set.md)
  handles continuous EEG data from EEGLAB much better
- Now using `whitening` package for whitening before SOBI ICA
- `select_epochs` for `eeg_ICA` objects fixed to correctly remove epochs
  from `signals`
- added tests for `filter.eeg_ICA` and `filter.eeg_tfr`
- fixed `filter.eeg_data` and `filter.eeg_evoked`
- `select_elecs` for `eeg_ICA` now correctly removes components from the
  unmixing matrix
- switched back to using `left_join` from `dplyr` in the `tag_events`
  function as an easy fix for sorting of events when tagging.
- fixed odd interaction between
  [`select()`](https://dplyr.tidyverse.org/reference/select.html) and
  [`validate_channels()`](https://craddm.github.io/eegUtils/reference/validate_channels.md)
  that reordered channel names in `chan_info`
- `eeg_decomp` now doing better job of filtering for `ssd` method
- various `tibble` related warnings and errors cleaned up.
- `method = "gam"` should now yield sensible results for
  [`geom_topo()`](https://craddm.github.io/eegUtils/reference/geom_topo.md)
- [`run_ICA()`](https://craddm.github.io/eegUtils/reference/run_ICA.md)
  and `eeg_decomp()` methods now return components ordered by percent
  variance explained (high to low)
- removed scaling of components in SOBI ICA method
- `browse_data().eeg_ica` grid res reduced
- `eeg_reference().eeg_epochs` was always average referencing, now
  fixed.
- cleaner code in
  [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  for biharmonic smooth
- [`compute_psd()`](https://craddm.github.io/eegUtils/reference/compute_psd.md)
  now demeans individual segments when doing Welch FFT; also no longer
  errors when only one segment per channel
- [`eeg_tfr()`](https://craddm.github.io/eegUtils/reference/eeg_tfr.md)
  internal structure modified to keep 4 dimensions even after averaging,
  for consistency
- `stat_summary_by_fill()` added to do averaging for raster plots
  effectively.
- [`convert_tfr()`](https://craddm.github.io/eegUtils/reference/convert_tfr.md)
  now properly returns converted data
- [`import_raw()`](https://craddm.github.io/eegUtils/reference/import_raw.md)
  fix for Brain Vision Analyzer files with date fields in the markers
- added `version` field to most objects

## eegUtils 0.5.0

#### Function changes

- Default settings for Infomax ICA changed to be similar to
  EEGLAB/Fieldtrip.
- Faster reading of bdf implemented. Old behaviour can be retained using
  `fast_bdf = FALSE` parameter to
  [`import_raw()`](https://craddm.github.io/eegUtils/reference/import_raw.md)
- `eeg_combine` now supports combining lists.
- `eeg_reference` now supports `eeg_epochs` and `eeg_evoked` objects.
  \### Internal changes / bug fixes
- `plot_butterfly` should now be faster again.
- Much faster reader for BDF implemented.
- `eeg_filter` added `demean` parameter so that removing channel/epoch
  means during filtering is now optional. Defaults to TRUE.
- Added artefact detection options for ICA objects
  - [`ar_acf()`](https://craddm.github.io/eegUtils/reference/ar_acf.md)
    checks for low autocorrelation
  - [`ar_chanfoc()`](https://craddm.github.io/eegUtils/reference/ar_chanfoc.md)
    checks for excessive channel focality (e.g. components that load
    mostly on one channel)
  - [`ar_trialfoc()`](https://craddm.github.io/eegUtils/reference/ar_trialfoc.md)
    checks for trial focality (components that load mostly on a few
    trials)
  - [`ar_eogcor()`](https://craddm.github.io/eegUtils/reference/ar_eogcor.md)
    checks for correlation with EOG channels
- `topoplot` plotting radius logic altered
- `compute_csd` now uses `eeg_reference` rather than `reref_eeg`
- Unmixing matrix for SSD decompositions fixed.
- `compute_tfr` reworked to be faster.
- Faster baseline correction implemented using Rcpp.
- Padding now used during `compute_tfr`, which greatly improves
  speed/accuracy; units may change but this is a change in scaling
  factor.
- `epoch_data` now uses a more robust way of determining time
  limits/samples to include in each epoch that no longer fails at some
  combinations of time limit and sampling rate
- `eeg_average` returns objects of class(`eeg_evoked`, `eeg_epochs`)
- Updated R requirement to \>= 3.2.0
- Updated rlang requirement to \>= 0.4.0
- `compute_psd` fix for single segment data
- updated use of `nest` and `unnest` in keeping with `tidyr 1.0.0`
- [`as.data.frame.eeg_tfr()`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_tfr.md)
  now fixed to output correctly

## eegUtils 0.4.0

#### Function changes

- Behaviour of `as.data.frame` methods changed.
  - `cond_label` parameter is deprecated
  - Information from the new `epochs` structure is now automatically
    added to the data.frame
- `tag_epochs` function added for labelling
- `run_ICA` now includes the `fICA` package version of `fastica`, and
  now supports running PCA before ICA.
- `apply_ICA` function added to remove ICA components.
- Behaviour of
  [`select_times()`](https://craddm.github.io/eegUtils/reference/select_times.md)
  changed to use exact supplied times rather than finding nearest times
  in the data.
- Some wrappers around `dplyr` functions added:
  - [`select()`](https://dplyr.tidyverse.org/reference/select.html) now
    works for selecting electrodes from `eeg_data` and `eeg_epochs`
    objects.
  - [`filter()`](https://dplyr.tidyverse.org/reference/filter.html)
    filters by time points or epochs from `eeg_data` and `eeg_epochs`
    objects.
  - [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) adds
    columns to the `signals` from `eeg_data` and `eeg_epochs` objects.
- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  now has a `groups` parameter that allows the possibility of facetting
  by event labels.
- more Biosemi montages added
- [`events()`](https://craddm.github.io/eegUtils/reference/events.md)
  function added to easily access and modify the events structure of all
  `eegUtils` objects.
- [`channels()`](https://craddm.github.io/eegUtils/reference/channels.md)
  function added to easily access and modify the chan_info structure of
  all `eegUtils` objects.
- [`epochs()`](https://craddm.github.io/eegUtils/reference/epochs.md)
  function added to access and modify epochs structure.
- [`ar_eogreg()`](https://craddm.github.io/eegUtils/reference/ar_eogreg.md)
  function added for removing eye movement activity using regression.
- [`eeg_filter()`](https://craddm.github.io/eegUtils/reference/eeg_filter.md)
  function added for a unified method of filtering using either FIR or
  IIR
  - [`eeg_filter()`](https://craddm.github.io/eegUtils/reference/eeg_filter.md)
    supports use of multiple threads/cores through the `future` package.
  - `iir_filt()` will be deprecated
- [`geom_topo()`](https://craddm.github.io/eegUtils/reference/geom_topo.md)
  extension for `ggplot2` added. Allows plotting of a topographical
  scalp maps using standard `ggplot2` functions.
- Default `grid_res` for topography related plots increased to 200.

#### Internal changes / bug fixes

- `data.table` now used in the following functions internally:
  - `reref_eeg()`
  - `iir_filt()`
  - [`eeg_FASTER()`](https://craddm.github.io/eegUtils/reference/ar_FASTER.md)
- `reref_eeg()` now correctly excludes electrodes as requested.
- `iir_filt()` now correctly respects epoch boundaries.
- New field `epochs` added to `eeg_data` and `eeg_epochs` objects.
- `chan_info` changes to make chan_info consistent across systems.
- Corrected scaling factor for PSD
- `eeg_combine` now checks and fixes `eeg_data` timing consistency
- `eeg_tfr` objects now use differently organised underlying matrices.
- `eeg_ICA` fixed unmixing matrices, which were transposed.

## eegUtils 0.3.0

#### Function changes

- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  now has a scaling parameter to scale the size of any lines or markers
  drawn on the plot.
- [`plot_tfr()`](https://craddm.github.io/eegUtils/reference/plot_tfr.md)
  function now useable, with baseline correction also added.
- [`rm_baseline()`](https://craddm.github.io/eegUtils/reference/rm_baseline.md)
  now handles `eeg_tfr` objects.
- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method
  added for `eeg_tfr` objects.
- [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md)
  function now available for use with Morlet wavelets.
- [`plot_psd()`](https://craddm.github.io/eegUtils/reference/plot_psd.md)
  now allows changing of FFT parameters (e.g. number of FFT points,
  segment length)
- Data selectors added for `eeg_tfr` objects
  (e.g. [`select_elecs()`](https://craddm.github.io/eegUtils/reference/select_elecs.md))

#### Internal changes/ bug fixes

- [`plot_timecourse()`](https://craddm.github.io/eegUtils/reference/plot_timecourse.md)
  overhauled to be S3 method.
- [`plot_butterfly()`](https://craddm.github.io/eegUtils/reference/plot_butterfly.md)
  reworked internally to be more efficient
- [`rm_baseline()`](https://craddm.github.io/eegUtils/reference/rm_baseline.md)
  simplified internally, reworked to use matrices; split to separate
  file.
- [`select_elecs()`](https://craddm.github.io/eegUtils/reference/select_elecs.md)
  now works for `eeg_evoked` objects
- `eeg_decomp` function in progress for performing SSD analyses
- Various methods added for TFR analyses
- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  improvements internally. Now offers potential for facetting.
- Some `dplyr` functions implemented internally for some objects.

## eegUtils 0.2.1

#### Function changes

- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  added highlights parameter to allow specific electrodes to be
  highlighted.
- [`run_ICA()`](https://craddm.github.io/eegUtils/reference/run_ICA.md)
  now offers extended Infomax and Fastica thanks to the `ica` package
- `plotly` is now a “suggested” package rather than a dependency
- [`plot_psd()`](https://craddm.github.io/eegUtils/reference/plot_psd.md)
  function added to calculate and plot the PSD for `eeg_epochs` and
  `eeg_data` objects
- [`plot_tfr()`](https://craddm.github.io/eegUtils/reference/plot_tfr.md)
  function added to handle `eeg_tfr` objects.
- [`erp_image()`](https://craddm.github.io/eegUtils/reference/erp_image.md)
  now works with `eeg_ICA` objects
- Generic print methods added for `eeg_epochs` and `eeg_data`
- [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md)
  function added to performed TFA on `eeg_epochs`
- [`epoch_data()`](https://craddm.github.io/eegUtils/reference/epoch_data.md)
  now warns if some events are not found rather than stops. Only stops
  if *no* events are found.

#### Internal changes/ bug fixes

- `reref_eeg()`
  - correctly excludes multiple named electrodes (i.e. passed as
    characters rather than numbers), where it previously silently
    failed.
  - no longer records the reference data in the `ref_data` field
- `tf_morlet` recoded to be called internally
- [`compute_psd()`](https://craddm.github.io/eegUtils/reference/compute_psd.md)
  - recoded to call
    [`welch_fft()`](https://craddm.github.io/eegUtils/reference/welch_fft.md)
    in order to support possibility of different FFT methods.
  - now drops the DC component (frequency 0)
- [`welch_fft()`](https://craddm.github.io/eegUtils/reference/welch_fft.md)
  internal function added
- [`eeg_downsample()`](https://craddm.github.io/eegUtils/reference/eeg_downsample.md)
  now makes sure epoch length is a multiple of the downsampling factor
  to avoid problems with timing jitter
- [`erp_image()`](https://craddm.github.io/eegUtils/reference/erp_image.md)
  is now an S3 method
- [`run_ICA()`](https://craddm.github.io/eegUtils/reference/run_ICA.md)
  - now returns source activations as a “signals” data frame, with
    component names
  - now returns correct unmixing matrix
- [`compute_csd()`](https://craddm.github.io/eegUtils/reference/compute_csd.md)
  - added.
  - computation of g-matrix and h-matrix refactored, spherical spline
    calculation altered accordingly
- [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md)
  - added
- [`eeg_FASTER()`](https://craddm.github.io/eegUtils/reference/ar_FASTER.md)
  now properly selects electrodes and epochs for removal

## eegUtils 0.2.0

#### Function changes

- [`as.data.frame.eeg_epochs()`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_epochs.md)
  now has a `cond_labels` parameter to select epochs with specific
  events and add the event label as an additional column.
- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  now drop the `sample` column.
- [`as.data.frame.eeg_ICA()`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_ICA.md)
  now has a `cond_labels` parameter to select epochs with specific
  events and add the event label as an additional column.
- `reref_eeg()` now removes reference channels from the data.
- [`eeg_FASTER()`](https://craddm.github.io/eegUtils/reference/ar_FASTER.md) -
  FASTER artefact rejection method now (mostly) implemented
  (*experimental*).

#### Internal changes / bug fixes

- [`plot_butterfly()`](https://craddm.github.io/eegUtils/reference/plot_butterfly.md)
  some `dplyr` use removed.
- `run_ica()` refactored SOBI method, JADE dependency removed.
- [`montage_check()`](https://craddm.github.io/eegUtils/reference/montage_check.md)
  command parses montage info when passed to
  [`electrode_locations()`](https://craddm.github.io/eegUtils/reference/electrode_locations.md)
- [`label_check()`](https://craddm.github.io/eegUtils/reference/label_check.md)
  added to help parse event labels
- [`proc_events()`](https://craddm.github.io/eegUtils/reference/proc_events.md)
  added to help parse event labels during
  [`select_epochs()`](https://craddm.github.io/eegUtils/reference/select_epochs.md)
  calls
- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  now tries to average/select across time/epochs before converting to
  long data, less memory use
- [`select_elecs()`](https://craddm.github.io/eegUtils/reference/select_elecs.md)
  also removes electrodes from chan_info
- [`select_epochs()`](https://craddm.github.io/eegUtils/reference/select_epochs.md)
  fixed bug where `events` and `timings` were inconsistent when using
  `keep = FALSE`
- Electrode/channel related functions (other than selection) now moved
  to `channel_management.r`
- New default electrode locations (347 locations in the 10-05 layout)
  provided

## eegUtils 0.1.15

#### Function changes

- [`eeg_evoked()`](https://craddm.github.io/eegUtils/reference/eeg_evoked.md)
  class introduced to hold ERPS
- [`eeg_ICA()`](https://craddm.github.io/eegUtils/reference/eeg_ICA.md)
  class introduced to hold ICA decompositions
- [`eeg_average()`](https://craddm.github.io/eegUtils/reference/eeg_average.md)
  function to calculate averages (e.g. ERPs) from `eeg_epochs` objects
- [`as.data.frame.eeg_evoked()`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_evoked.md)
  introduced to handle conversion of eeg_evoked objects to data frames.

#### Internal changes / bug fixes

- [`compute_psd()`](https://craddm.github.io/eegUtils/reference/compute_psd.md)
  function development, converted to S3method.
- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  properly checks for existing chan_info in `eeg_data` objects
- [`plot_timecourse()`](https://craddm.github.io/eegUtils/reference/plot_timecourse.md)
  and
  [`plot_butterfly()`](https://craddm.github.io/eegUtils/reference/plot_butterfly.md)
  modified to deal with `eeg_evoked` objects.
- [`plot_butterfly()`](https://craddm.github.io/eegUtils/reference/plot_butterfly.md)
  updated to better handle data frames
- [`topoplot.eeg_ICA()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  added to make topolots from ICA components
- [`rm_baseline()`](https://craddm.github.io/eegUtils/reference/rm_baseline.md)
  reworked as S3 method and to be faster and much less memory intensive.
- [`plot_butterfly()`](https://craddm.github.io/eegUtils/reference/plot_butterfly.md)
  converted to S3 method.
- Initial commits for addition of Morlet wavelet time-frequency analysis
- Initial commits for statisical comparisons added

## eegUtils 0.1.14

#### Function changes

- [`eeg_downsample()`](https://craddm.github.io/eegUtils/reference/eeg_downsample.md)
  function added to downsample EEG data by an integer factor.
- [`tag_events()`](https://craddm.github.io/eegUtils/reference/tag_events.md)
  function added to give labels to event codes.
- [`list_events()`](https://craddm.github.io/eegUtils/reference/list_events.md)
  added to display unique event codes and their associated labels.
- [`select_epochs()`](https://craddm.github.io/eegUtils/reference/select_epochs.md)
  now allows selection of epochs by event code or event label.
- [`erp_raster()`](https://craddm.github.io/eegUtils/reference/erp_raster.md) -
  plot ERPs across the scalp as an ERP image
- [`eeg_combine()`](https://craddm.github.io/eegUtils/reference/eeg_combine.md) -
  combine multiple `eeg_data` or `eeg_epochs` objects into one

#### Internal changes/ bug fixes

- [`eeg_epochs()`](https://craddm.github.io/eegUtils/reference/eeg_epochs.md)
  now also handles downsampled data appropriately.
- [`select_times()`](https://craddm.github.io/eegUtils/reference/select_times.md)
  no longer leaves “epoch” column in `eeg_epochs` objects.
- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  now calls a separate function (`gam_topo()`) to create GAM smooths
- [`browse_data()`](https://craddm.github.io/eegUtils/reference/browse_data.md)
  major speed-ups, no longer converts to long format until necessary.
  Converted to S3method.
- [`interactive_scalp()`](https://craddm.github.io/eegUtils/reference/interactive_scalp.md)
  fixed plotting of individual electrodes

## eegUtils 0.1.13

#### Function changes

- [`interp_elecs()`](https://craddm.github.io/eegUtils/reference/interp_elecs.md)
  function to perform spherical spline interpolation of individual
  electrodes
- `eeg_ar_thresh()` simple absolute value thresholding added
- [`plot_electrodes()`](https://craddm.github.io/eegUtils/reference/plot_electrodes.md)
  Produces a 2D or interactive 3D plot of electrode locations

## eegUtils 0.1.12

#### Function changes

- `iir_filt()` now also filters reference channels
- `load_set()` command added to load EEGLAB .set files

#### Internal changes

- Converted
  [`select_times()`](https://craddm.github.io/eegUtils/reference/select_times.md)
  to an S3 generic method
  - `select_times.eeg_data`
  - `select_times.eeg_epochs`
- Converted `iir_filt()` to an S3 generic method
  - `iir_filt.eeg_data`
  - `iir_filt.eeg_epochs`

## eegUtils 0.1.11

- Added a `NEWS.md` file to track changes to the package.

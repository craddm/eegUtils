# eegUtils (development version)

### Function changes
- `plot_timecourse()` now requires facets and mappings to be explicitly stated during the call, rather than added afterwards. This allows it to use weighted averages when required.
- `plot_timecourse.eeg_tfr()` now defaults to no baseline correction, and takes a `type` argument for baseline type to be specified.
- Added new `plot_gfp()` function for calculating and plotting Global Field Power.
- `topoplot()` now supports plotting of multiple timepoints - pass a list of times to the `time_lim` argument.
- `topoplot()` now supports custom fill titles through the `fill_title` argument, and automatically switches where necessary (e.g. now says "Power mV^2"). Fill titles are now also centred.
- `browse_data.eeg_ICA()` now provides the option to select components for rejection, and returns a character vector of selected components.
- `view_ica()` now allows you to select components for rejection, to double-click on topographies to inspect them individually, and to return cleaned data.
- `eeg_average()` now supports averaging over conditions in `eeg_evoked` files
- `eeg_average()` now records weights - the number of epochs that went into an average - and uses those in subsequent steps where possible for `eeg_epochs`/`eeg_evoked` objects. Currently only partially implemented for `eeg_tfr` objects.
- `compute_tfr()` now has an argument `trim_edges` which allows users to switch off automatic removal of epoch edges after transformation. Defaults to TRUE.
- `rm_baseline()` now adds a record of the baseline period 

### Internal changes/bug fixes
- Recoded `faster_epochs()` to no longer use `data.table`.
- added `parse_cycles()` function for use during `compute_tfr()`
- Using `electrode_locations()` on a `data.frame` would return a data frame with the electrode names in full upper case. Now returns with the electrodes in their original case.
- `eeg_evoked` objects should now contain reference information.
- `ar_for_ica` file started, moving some functions from `artefact_rejection.R`
- `select_times` prevented from converting single column data.frames to vectors after subsetting.
- `plot_timecourse.eeg_tfr` now correctly uses the `freq_range` parameter.
- Improvements to internal processing logic in `view_ica()` to improve performance.
- `import_raw(..., fast_bdf = TRUE)` will now discard Annotations rather than fail to import BDF files with Annotations.
- removed dependency on `Matrix` - now using `qr()$rank` in `run_ICA` rather than using `Matrix::rankMatrix()` to determine rank of input signals.
- Added numerous `dropped_aes` variables to the custom `ggplot2` `stat_` functions for compatibility with `ggplot2 3.4.0`.
- Replaced `size` aesthetic with `linewidth` for compatibility with `ggplot2 3.4.0`.

# eegUtils 0.7.0

### Function changes
- Add `imax` method to `run_ICA()`. This allows use of the `infomax` ICA algorithm from the `infomax` package, which is a reimplementation of the extended-Infomax algorithm used in the `EEGLAB` Matlab toolbox.
- `erp_scalp()` and `interactive_scalp()` should now appropriately use channel locations included in the data.
- More informative messages when using `compute_tfr()`.
- `plot_tfr()` now applies baseline correction on a single-trial basis where possible, which may show different results when using non-linear baseline correction (e.g. `divide` or `dB`)
- `topoplot()` now allows you to provide multiple component numbers when plotting from an `eeg_ICA` object, and will automatically produce an appropriately facetted plot. 
- `topoplot()` now has a `k` parameter to control the smoothing when using `method = "gam"`.
- added additional `demo_spatial` data from a spatial cueing experiment.
- `plot_difference()` function added for plotting ERP difference waves. Only currently handles two levels.
- Added `hanning` taper support for `compute_tfr`. Note that the scaling factors used for all `compute_tfr` calculations have been adjusted, so the exact numerical values returned will change. However, this is just a scaling factor - the relative distances between values remained unchanged.
- `ar_FASTER()` has experimental support for `eeg_group` objects when those objects are `eeg_evoked` groups. It does not perform rejection but reports how many times each participants data breaks a threshold for a number of measures.
- `import_raw()` default `participant_id` is now changed to `NA` instead of `character(1)`, to promote better use with `eeg_combine()`.
- `eeg_combine()` will now refuse to combine objects where `participant_id` is missing (i.e. is `NA`), and warn when combining objects with the previous default value "". This is to prevent accidentally treating data from different participants as being from the same participant.
- `topoplot()` now provides informative messages about the head radius used for plotting. The default calculation of `r` when using `interp_limit = "head"` has changed and should now be set at the outermost electrode's position + a 10% buffer. Smaller `r` can be set manually.
- added `get_participant_id`, `set_participant_id`, `get_recording` and `set_recording` to interact with the `epochs` metadata in each `eegUtils` object.
- `eeg_combine()` handles `eeg_tfr()` objects better, now returns an error when trying to combine single-trial data.
- `plot_timecourse()` now handles `eeg_tfr()` objects.

### Internal changes / bug fixes

- `erp_scalp()` and `interactive_scalp()` now use cleaner evaluation of the `colour` argument using `rlang` 
- Fixed adding CSD as new reference when using `compute_csd()`
- When using a scaling number of cycles in `compute_tfr()`, they will now also use the `spacing` parameter to determine `log` or `linear` scaling.
- Fixed bug with `eeg_average()` used on `c("eeg_group", "eeg_tfr")` objects.
- Fixed bug with incorrect number of epochs calculated when epoching `eeg_data` objects if there were multiple "target" triggers appearing in an epoch.
- Fixed error with `erp_raster()` when `anat_order == FALSE`
- new `epoch_queries` file for functions for setting and getting the `epochs` structure
- new `check_items` file for functions that check for consistency of various structures
- `eeg_combine()` with `eeg_tfr()` no longer drops single dimensions, which was causing issues when there was only one channel or epoch in the object.
- `plot_timecourse.eeg_tfr` now correctly passes baseline period to `rm_baseline()`
- `filter.eeg_tfr()` was sometimes dropping single dimensions when using `abind::asub()`, now fixed

# eegUtils 0.6.3

### Function changes
- Added log spaced frequencies to `compute_tfr()`, with the new `spacing` argument. `plot_tfr` automatically detects the spacing and plots the figure appropriately.
- Added `na.rm` option to `erp_image()` to either keep or plot trials with NA values due to smoothing. By default they'll be removed.
- Added more informative messages for `compute_psd()`.

### Internal changes / bug fixes

- Some minor documentation fixes.
- `plot_tfr()` error when selecting a specific frequency range fixed.
- Switched to new style of `vdiffr` tests
- fixed `erp_image()` smoothing over time instead of epochs.

# eegUtils 0.6.2

### Function changes
- added support for `EEGLAB` .set files saved in newer Matlab file formats.
- changed first argument of `eeg_filter()` to `data` instead of `.data`
- added some more informative user messages for importing data and adding electrode locations.
- `erp_image()` now supports `eeg_tfr` objects.

### Internal changes / bug fixes

- When combining three or more continuous `eeg_data` objects, `eeg_combine()` would substantially undercorrect the timing of events in the third file - this is now fixed.
- `groups` parameter for `topoplot()` now correctly passed for all types of object.
- `stat_tests.R` file removed, will be reimplemented elsewhere
- Long standing issues with import of channel locations from EEGLAB files hopefully fixed...
- `rm_baseline()` for `eeg_evoked` no longer uses `data.table`
- `as.data.frame.eeg_evoked()` handles grouped data better.
- `import_set()` now handles all EEGLAB formats better.


# eegUtils 0.6.1

### Function changes

- `ar_eogcor` now has a `bipolarize` argument which can be set to false when the HEOG/VEOG channels are already bipolarized.
- added some new `ggplot2` based functions for topoplotting and adding contours   
    - `stat_scalpcontours()`
- `geom_topo()` now has contours
- added errors with when attempting to use `compute_psd` or `compute_tfr` on `eeg_group` objects.
- `compute_tfr()` now works better with `eeg_evoked` objects that contain multiple conditions.


### Internal changes / bug fixes

- `eeg_reference` now handles multiple reference channels better on rereferencing
- `get_scalpmap` handles `eeg_ICA` components better when there are channels with no locations
- Travis-CI removed.
- `cart_to_spherical` coord flipping bug fixed (hopefully...)
- `filter` now converts to tibble internally and does not coerce `signals` to a vector when there is only one channel.
- added copyright info to `summary_contour` file
- `eeg_combine.eeg_evoked` made to behave more consistently when creating grouped data
- all `as.data.frame()` functions moved to `df_converters.r`

# eegUtils 0.6.0

IMPORTANT: There have been some changes to the logic of the `topoplot()` that may make their appearance quite different. Specifically, these changes are to the way the underlying interpolation grid is calculated and to how things like the diameter of the cartoon head is calculated. These changes often lead to different minimum or maximum amplitudes across the image, and thus changes in the appearance of the plot due to different scales- don't be alarmed!

### Function changes
- `epoch_data` baseline correction now defaults to *no* baseline correction
- Added `filter` method for `eeg_tfr` objects
- `fit_glm()` overhauled. Now far faster and allows specification of models using standard R formulae.
- New `eeg_lm` class introduced for output of `fit_glm()`.
- `plot_butterfly.eeg_lm()` method added
- `as.data.frame()` methods have been added for `eeg_lm` objects.
- `view_ica()` Shiny viewer for `eeg_ICA` and `eeg_decomp` objects added.
- `view_artefacts()`Shiny viewer for channel and epoch stats added.
- `plot_timecourse()` now takes a mapping argument, which allows use of `ggplot2` `aes()` mappings
- `eeg_average.eeg_tfr()` now follows behaviour of other `eeg_average()` methods in respecting the `epochs` structure.  
- `compute_itc()` added for computing inter-trial coherence from `eeg_tfr` objects.
- `cols` added to `eeg_average.eeg_tfr`
- `eeg_combine.tfr_average()` added to handle pre-averaged `eeg_tfr` objects
- `compute_tfr()` now allows non-constant number of cycles
- `compute_tfr()` now uses a different scaling factor, so raw units should now be microvolts-squared. 
- added `import_erplab()` function
- `plot_timecourse()` now allows CIs for `eeg_group` objects.

### Internal changes / bug fixes
- `plot_tfr()` now always drops NA/NaN values and averages appropriately over electrodes and conditions.
- `import_set()` handles continuous EEG data from EEGLAB much better
- Now using `whitening` package for whitening before SOBI ICA
- `select_epochs` for `eeg_ICA` objects fixed to correctly remove epochs from `signals`
- added tests for `filter.eeg_ICA` and `filter.eeg_tfr`
- fixed `filter.eeg_data` and `filter.eeg_evoked`
- `select_elecs` for `eeg_ICA` now correctly removes components from the unmixing matrix
- switched back to using `left_join` from `dplyr` in the `tag_events` function as an easy fix for sorting of events when tagging.
- fixed odd interaction between `select()` and `validate_channels()` that reordered channel names in `chan_info`
- `eeg_decomp` now doing better job of filtering for `ssd` method
- various `tibble` related warnings and errors cleaned up.
- `method = "gam"` should now yield sensible results for `geom_topo()`
- `run_ICA()` and `eeg_decomp()` methods now return components ordered by percent variance explained (high to low)
- removed scaling of components in SOBI ICA method
- `browse_data().eeg_ica` grid res reduced
- `eeg_reference().eeg_epochs` was always average referencing, now fixed.
- cleaner code in `topoplot()` for biharmonic smooth
- `compute_psd()` now demeans individual segments when doing Welch FFT; also no longer errors when only one segment per channel
- `eeg_tfr()` internal structure modified to keep 4 dimensions even after averaging, for consistency
- `stat_summary_by_fill()` added to do averaging for raster plots effectively.
- `convert_tfr()` now properly returns converted data
- `import_raw()` fix for Brain Vision Analyzer files with date fields in the markers
- added `version` field to most objects


# eegUtils 0.5.0

### Function changes
- Default settings for Infomax ICA changed to be similar to EEGLAB/Fieldtrip.
- Faster reading of bdf implemented. Old behaviour can be retained using `fast_bdf = FALSE` parameter to `import_raw()`
- `eeg_combine` now supports combining lists.
- `eeg_reference` now supports `eeg_epochs` and `eeg_evoked` objects.
### Internal changes / bug fixes
- `plot_butterfly` should now be faster again.
- Much faster reader for BDF implemented.
- `eeg_filter` added `demean` parameter so that removing channel/epoch means during filtering is now optional. Defaults to TRUE.
- Added artefact detection options for ICA objects
    - `ar_acf()` checks for low autocorrelation
    - `ar_chanfoc()` checks for excessive channel focality (e.g. components that load mostly on one channel)
    - `ar_trialfoc()` checks for trial focality (components that load mostly on a few trials)
    - `ar_eogcor()` checks for correlation with EOG channels
- `topoplot` plotting radius logic altered 
- `compute_csd` now uses `eeg_reference` rather than `reref_eeg`
- Unmixing matrix for SSD decompositions fixed.
- `compute_tfr` reworked to be faster.
- Faster baseline correction implemented using Rcpp.
- Padding now used during `compute_tfr`, which greatly improves speed/accuracy; units may change but this is a change in scaling factor.
- `epoch_data` now uses a more robust way of determining time limits/samples to include in each epoch that no longer fails at some combinations of time limit and sampling rate
- `eeg_average` returns objects of class(`eeg_evoked`, `eeg_epochs`)
- Updated R requirement to >= 3.2.0
- Updated rlang requirement to >= 0.4.0
- `compute_psd` fix for single segment data
- updated use of `nest` and `unnest` in keeping with `tidyr 1.0.0`
- `as.data.frame.eeg_tfr()` now fixed to output correctly

# eegUtils 0.4.0

### Function changes
- Behaviour of `as.data.frame` methods changed.
    - `cond_label` parameter is deprecated
    - Information from the new `epochs` structure is now automatically added to the data.frame
- `tag_epochs` function added for labelling
- `run_ICA` now includes the `fICA` package version of `fastica`, and now supports running PCA before ICA.
- `apply_ICA` function added to remove ICA components.
- Behaviour of `select_times()` changed to use exact supplied times rather than finding nearest times in the data.
- Some wrappers around `dplyr` functions added:
    - `select()` now works for selecting electrodes from `eeg_data` and `eeg_epochs` objects.
    - `filter()` filters by time points or epochs from `eeg_data` and `eeg_epochs` objects.
    - `mutate()` adds columns to the `signals` from `eeg_data` and `eeg_epochs` objects.
- `topoplot()` now has a `groups` parameter that allows the possibility of facetting by event labels.
- more Biosemi montages added
- `events()` function added to easily access and modify the events structure of all `eegUtils` objects.
- `channels()` function added to easily access and modify the chan_info structure of all `eegUtils` objects.
- `epochs()` function added to access and modify epochs structure.
- `ar_eogreg()` function added for removing eye movement activity using regression.
- `eeg_filter()` function added for a unified method of filtering using either FIR or IIR
    - `eeg_filter()` supports use of multiple threads/cores through the `future` package.
    - `iir_filt()` will be deprecated
- `geom_topo()` extension for `ggplot2` added. Allows plotting of a topographical scalp maps using standard `ggplot2` functions.
- Default `grid_res` for topography related plots increased to 200.

### Internal changes / bug fixes
- `data.table` now used in the following functions internally:
    - `reref_eeg()`
    - `iir_filt()` 
    - `eeg_FASTER()`
- `reref_eeg()` now correctly excludes electrodes as requested.
- `iir_filt()` now correctly respects epoch boundaries.
- New field `epochs` added to `eeg_data` and `eeg_epochs` objects.
- `chan_info` changes to make chan_info consistent across systems. 
- Corrected scaling factor for PSD
- `eeg_combine` now checks and fixes `eeg_data` timing consistency
- `eeg_tfr` objects now use differently organised underlying matrices.
- `eeg_ICA` fixed unmixing matrices, which were transposed. 

# eegUtils 0.3.0

### Function changes
- `topoplot()` now has a scaling parameter to scale the size of any lines or markers drawn on the plot.
- `plot_tfr()` function now useable, with baseline correction also added.
- `rm_baseline()` now handles `eeg_tfr` objects.
- `as.data.frame()` method added for `eeg_tfr` objects.
- `compute_tfr()` function now available for use with Morlet wavelets.
- `plot_psd()` now allows changing of FFT parameters (e.g. number of FFT points, segment length)
- Data selectors added for `eeg_tfr` objects (e.g. `select_elecs()`)

### Internal changes/ bug fixes
- `plot_timecourse()` overhauled to be S3 method.
- `plot_butterfly()` reworked internally to be more efficient
- `rm_baseline()` simplified internally, reworked to use matrices; split to separate file.
- `select_elecs()` now works for `eeg_evoked` objects
- `eeg_decomp` function in progress for performing SSD analyses
- Various methods added for TFR analyses
- `topoplot()` improvements internally. Now offers potential for facetting.
- Some `dplyr` functions implemented internally for some objects.

# eegUtils 0.2.1

### Function changes
- `topoplot()` added highlights parameter to allow specific electrodes to be highlighted.
- `run_ICA()` now offers extended Infomax and Fastica thanks to the `ica` package
- `plotly` is now a "suggested" package rather than a dependency
- `plot_psd()` function added to calculate and plot the PSD for `eeg_epochs` and `eeg_data` objects
- `plot_tfr()` function added to handle `eeg_tfr` objects.
- `erp_image()` now works with `eeg_ICA` objects
- Generic print methods added for `eeg_epochs` and `eeg_data`
- `compute_tfr()` function added to performed TFA on `eeg_epochs`
- `epoch_data()` now warns if some events are not found rather than stops. Only stops if *no* events are found.

### Internal changes/ bug fixes
- `reref_eeg()` 
    - correctly excludes multiple named electrodes (i.e. passed as characters rather than numbers), where it previously silently failed.
    - no longer records the reference data in the `ref_data` field
- `tf_morlet` recoded to be called internally
- `compute_psd()` 
    - recoded to call `welch_fft()` in order to support possibility of different FFT methods.
    - now drops the DC component (frequency 0)
- `welch_fft()` internal function added
- `eeg_downsample()` now makes sure epoch length is a multiple of the downsampling factor to avoid problems with timing jitter
- `erp_image()` is now an S3 method
- `run_ICA()` 
    - now returns source activations as a "signals" data frame, with component names
    - now returns correct unmixing matrix
- `compute_csd()` 
  - added.
  - computation of g-matrix and h-matrix refactored, spherical spline calculation altered accordingly
- `compute_tfr()` 
  - added
- `eeg_FASTER()` now properly selects electrodes and epochs for removal

# eegUtils 0.2.0

### Function changes
- `as.data.frame.eeg_epochs()` now has a `cond_labels` parameter to select epochs with specific events and add the event label as an additional column.
- `as.data.frame()` methods now drop the `sample` column.
- `as.data.frame.eeg_ICA()` now has a `cond_labels` parameter to select epochs with specific events and add the event label as an additional column.
- `reref_eeg()` now removes reference channels from the data.
- `eeg_FASTER()` - FASTER artefact rejection method now (mostly) implemented (*experimental*).

### Internal changes / bug fixes
- `plot_butterfly()` some `dplyr` use removed.
- `run_ica()` refactored SOBI method, JADE dependency removed.
- `montage_check()` command parses montage info when passed to `electrode_locations()`
- `label_check()` added to help parse event labels
- `proc_events()` added to help parse event labels during `select_epochs()` calls
- `topoplot()` now tries to average/select across time/epochs before converting to long data, less memory use
- `select_elecs()` also removes electrodes from chan_info
- `select_epochs()` fixed bug where `events` and `timings` were inconsistent when using `keep = FALSE`
- Electrode/channel related functions (other than selection) now moved to `channel_management.r`
- New default electrode locations (347 locations in the 10-05 layout) provided

# eegUtils 0.1.15

### Function changes
- `eeg_evoked()` class introduced to hold ERPS
- `eeg_ICA()` class introduced to hold ICA decompositions
- `eeg_average()` function to calculate averages (e.g. ERPs) from `eeg_epochs` objects
- `as.data.frame.eeg_evoked()` introduced to handle conversion of eeg_evoked objects to data frames.

### Internal changes / bug fixes
- `compute_psd()` function development, converted to S3method.
- `topoplot()` properly checks for existing chan_info in `eeg_data` objects
- `plot_timecourse()` and `plot_butterfly()` modified to deal with `eeg_evoked` objects.
- `plot_butterfly()` updated to better handle data frames
- `topoplot.eeg_ICA()` added to make topolots from ICA components
- `rm_baseline()` reworked as S3 method and to be faster and much less memory intensive.
- `plot_butterfly()` converted to S3 method.
- Initial commits for addition of Morlet wavelet time-frequency analysis
- Initial commits for statisical comparisons added

# eegUtils 0.1.14

### Function changes
- `eeg_downsample()` function added to downsample EEG data by an integer factor.
- `tag_events()` function added to give labels to event codes.
- `list_events()` added to display unique event codes and their associated labels.
- `select_epochs()` now allows selection of epochs by event code or event label.
- `erp_raster()` - plot ERPs across the scalp as an ERP image
- `eeg_combine()` - combine multiple `eeg_data` or `eeg_epochs` objects into one

### Internal changes/ bug fixes
- `eeg_epochs()` now also handles downsampled data appropriately.
- `select_times()` no longer leaves "epoch" column in `eeg_epochs` objects.
- `topoplot()` now calls a separate function (`gam_topo()`) to create GAM smooths
- `browse_data()` major speed-ups, no longer converts to long format until necessary. Converted to S3method.
- `interactive_scalp()` fixed plotting of individual electrodes

# eegUtils 0.1.13

### Function changes
 - `interp_elecs()` function to perform spherical spline interpolation of individual electrodes
 - `eeg_ar_thresh()` simple absolute value thresholding added
 - `plot_electrodes()` Produces a 2D or interactive 3D plot of electrode locations

# eegUtils 0.1.12

### Function changes
 - `iir_filt()` now also filters reference channels
 - `load_set()` command added to load EEGLAB .set files

### Internal changes

* Converted `select_times()` to an S3 generic method
  * `select_times.eeg_data`
  * `select_times.eeg_epochs`

* Converted `iir_filt()` to an S3 generic method
  * `iir_filt.eeg_data`
  * `iir_filt.eeg_epochs`

# eegUtils 0.1.11

* Added a `NEWS.md` file to track changes to the package.

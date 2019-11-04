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
- `epoch_data` now uses a more robust way of determing time limits/samples to include in each epoch that no longer fails at some combinations of time limit and sampling rate
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

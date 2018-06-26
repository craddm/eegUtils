# eegUtils 0.2.0.9000

### Function changes

### Internal changes/ bug fixes
- `reref_eeg()` correctly excludes multiple named electrodes (i.e. passed as characters rather than numbers), where it previously silently failed.
- `tf_morlet` recoded to be called internally
- `compute_psd` recoded to call `welch_fft()` in order to support possibility of different FFT methods.
- `welch_fft()` internal function added

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

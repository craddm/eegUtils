# eegUtils 0.1.14

### Function changes
- `eeg_downsample()` function added to downsample EEG data by an integer factor.
- `tag_events()` function added to give labels to event codes.
- `list_events()` added to display unique event codes and their associated labels.
- `select_epochs()` now allows selection of epochs by event code or event label.

### Internal changes
- `eeg_epochs()` now also handles downsampled data appropriately.

# eegUtils 0.1.13

### Function changes
 - `interp_elecs()` function to perform spherical spline interpolation of individual electrodes.
 - `eeg_ar_thresh()` simple absolute value thresholding added.
 - `plot_electrodes()` Produces a 2D or interactive 3D plot of electrode locations.

# eegUtils 0.1.12

### Function changes
 - `iir_filt()` now also filters reference channels.
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

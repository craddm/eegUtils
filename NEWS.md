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

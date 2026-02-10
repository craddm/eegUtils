# Package index

## IO

Functions for importing data or electrode information into R.

- [`import_chans()`](https://craddm.github.io/eegUtils/reference/import_chans.md)
  : Import channel locations from various file formats

- [`import_erplab()`](https://craddm.github.io/eegUtils/reference/import_erplab.md)
  : Import from ERPLAB .erp files

- [`import_ft()`](https://craddm.github.io/eegUtils/reference/import_ft.md)
  : Import Fieldtrip files

- [`import_raw()`](https://craddm.github.io/eegUtils/reference/import_raw.md)
  : Function for reading raw data.

- [`import_set()`](https://craddm.github.io/eegUtils/reference/import_set.md)
  :

  Load `EEGLAB` .set files

- [`export_bva()`](https://craddm.github.io/eegUtils/reference/export_bva.md)
  : Export continuous data in Brain Vision Analyzer format

## Processing

Functions for pre-processing and processing data

- [`apply_ica()`](https://craddm.github.io/eegUtils/reference/apply_ica.md)
  : Recreate channel timecourses from ICA decompositions.

- [`compute_csd()`](https://craddm.github.io/eegUtils/reference/compute_csd.md)
  : Convert to Current Source Density

- [`eeg_average()`](https://craddm.github.io/eegUtils/reference/eeg_average.md)
  : Calculate averages (e.g. event-related potentials) for single
  datasets

- [`eeg_combine()`](https://craddm.github.io/eegUtils/reference/eeg_combine.md)
  :

  Combine `eegUtils` objects

- [`eeg_decompose()`](https://craddm.github.io/eegUtils/reference/eeg_decompose.md)
  : Generalized eigenvalue decomposition based methods for EEG data

- [`eeg_downsample()`](https://craddm.github.io/eegUtils/reference/eeg_downsample.md)
  : Downsampling EEG data

- [`eeg_filter()`](https://craddm.github.io/eegUtils/reference/eeg_filter.md)
  : Filter EEG data

- [`eeg_reference()`](https://craddm.github.io/eegUtils/reference/eeg_reference.md)
  : Referencing

- [`electrode_locations()`](https://craddm.github.io/eegUtils/reference/electrode_locations.md)
  : Get standard electrode locations

- [`epoch_data()`](https://craddm.github.io/eegUtils/reference/epoch_data.md)
  : Create epochs from EEG data

- [`interp_elecs()`](https://craddm.github.io/eegUtils/reference/interp_elecs.md)
  : Channel interpolation

- [`rm_baseline()`](https://craddm.github.io/eegUtils/reference/rm_baseline.md)
  : Baseline correction

- [`run_ICA()`](https://craddm.github.io/eegUtils/reference/run_ICA.md)
  : Independent Component Analysis for EEG data

- [`tag_events()`](https://craddm.github.io/eegUtils/reference/tag_events.md)
  : Tag events

## Artefact rejection

Functions for artefact rejection.

- [`ar_FASTER()`](https://craddm.github.io/eegUtils/reference/ar_FASTER.md)
  [`eeg_FASTER()`](https://craddm.github.io/eegUtils/reference/ar_FASTER.md)
  : FASTER EEG artefact rejection
- [`ar_acf()`](https://craddm.github.io/eegUtils/reference/ar_acf.md) :
  Detect low autocorrelation of ICA components
- [`ar_chanfoc()`](https://craddm.github.io/eegUtils/reference/ar_chanfoc.md)
  : Detect high channel focality of ICA components
- [`ar_eogcor()`](https://craddm.github.io/eegUtils/reference/ar_eogcor.md)
  : Detect high component correlation with eye channels
- [`ar_eogreg()`](https://craddm.github.io/eegUtils/reference/ar_eogreg.md)
  : Remove EOG using regression
- [`ar_thresh()`](https://craddm.github.io/eegUtils/reference/ar_thresh.md)
  : Simple absolute value thresholding
- [`ar_trialfoc()`](https://craddm.github.io/eegUtils/reference/ar_trialfoc.md)
  : Detect high trial focality of ICA components
- [`epoch_stats()`](https://craddm.github.io/eegUtils/reference/epoch_stats.md)
  : Epoch statistics
- [`channel_stats()`](https://craddm.github.io/eegUtils/reference/channel_stats.md)
  : Channel statistics
- [`view_artefacts()`](https://craddm.github.io/eegUtils/reference/view_artefacts.md)
  : Artefact browser

## Selection

Functions for selecting subsets of data

- [`select_elecs()`](https://craddm.github.io/eegUtils/reference/select_elecs.md)
  : Select electrodes from a given dataset
- [`select_epochs()`](https://craddm.github.io/eegUtils/reference/select_epochs.md)
  : Select epochs
- [`select_freqs()`](https://craddm.github.io/eegUtils/reference/select_freqs.md)
  : Select frequencies
- [`select_times()`](https://craddm.github.io/eegUtils/reference/select_times.md)
  : Select timerange
- [`reexports`](https://craddm.github.io/eegUtils/reference/reexports.md)
  [`filter`](https://craddm.github.io/eegUtils/reference/reexports.md)
  [`select`](https://craddm.github.io/eegUtils/reference/reexports.md)
  [`mutate`](https://craddm.github.io/eegUtils/reference/reexports.md)
  [`rename`](https://craddm.github.io/eegUtils/reference/reexports.md)
  [`fortify`](https://craddm.github.io/eegUtils/reference/reexports.md)
  : Objects exported from other packages

## Plotting

Functions for plotting data

- [`browse_data()`](https://craddm.github.io/eegUtils/reference/browse_data.md)
  : Browse EEG data
- [`erp_image()`](https://craddm.github.io/eegUtils/reference/erp_image.md)
  : Plot ERP images
- [`erp_raster()`](https://craddm.github.io/eegUtils/reference/erp_raster.md)
  : ERP raster plot
- [`erp_scalp()`](https://craddm.github.io/eegUtils/reference/erp_scalp.md)
  : Plot event-related potentials using a scalp based layout
- [`geom_topo()`](https://craddm.github.io/eegUtils/reference/geom_topo.md)
  : Create a topographical plot
- [`get_scalpmap()`](https://craddm.github.io/eegUtils/reference/get_scalpmap.md)
  : Calculate an interpolated scalpmap
- [`interactive_scalp()`](https://craddm.github.io/eegUtils/reference/interactive_scalp.md)
  : Interactive scalp maps
- [`plot_butterfly()`](https://craddm.github.io/eegUtils/reference/plot_butterfly.md)
  : Create a butterfly plot from timecourse data
- [`plot_difference()`](https://craddm.github.io/eegUtils/reference/plot_difference.md)
  : Plot ERP difference waves
- [`plot_electrodes()`](https://craddm.github.io/eegUtils/reference/plot_electrodes.md)
  : Plot electrode locations
- [`plot_gfp()`](https://craddm.github.io/eegUtils/reference/plot_gfp.md)
  : Plot Global Field Power of EEG Signals
- [`plot_psd()`](https://craddm.github.io/eegUtils/reference/plot_psd.md)
  : Plot Power Spectral Density
- [`plot_tfr()`](https://craddm.github.io/eegUtils/reference/plot_tfr.md)
  : Time-frequency plot
- [`plot_timecourse()`](https://craddm.github.io/eegUtils/reference/plot_timecourse.md)
  : Plot one-dimensional timecourse data.
- [`stat_scalpmap()`](https://craddm.github.io/eegUtils/reference/stat_scalpmap.md)
  [`geom_head()`](https://craddm.github.io/eegUtils/reference/stat_scalpmap.md)
  [`geom_mask()`](https://craddm.github.io/eegUtils/reference/stat_scalpmap.md)
  [`geom_ears()`](https://craddm.github.io/eegUtils/reference/stat_scalpmap.md)
  [`geom_channels()`](https://craddm.github.io/eegUtils/reference/stat_scalpmap.md)
  : Create an interpolated scalp surface
- [`stat_scalpcontours()`](https://craddm.github.io/eegUtils/reference/stat_scalpcontours.md)
  : Create an interpolated scalp surface
- [`topoplot()`](https://craddm.github.io/eegUtils/reference/topoplot.md)
  : Topographical Plotting Function for EEG
- [`view_ica()`](https://craddm.github.io/eegUtils/reference/view_ica.md)
  : EEG decomposition viewer

## Frequency analysis

Functions related to (time-)frequency analysis

- [`compute_itc()`](https://craddm.github.io/eegUtils/reference/compute_ITC.md)
  : Calculate inter-trial coherence
- [`compute_psd()`](https://craddm.github.io/eegUtils/reference/compute_psd.md)
  : Compute power spectral density
- [`compute_tfr()`](https://craddm.github.io/eegUtils/reference/compute_tfr.md)
  : Compute Time-Frequency representation of EEG data

## Converters

Functions for converting objects to data.frames

- [`as.data.frame(`*`<eeg_ICA>`*`)`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_ICA.md)
  :

  Convert `eeg_ICA` object to data frame

- [`as.data.frame(`*`<eeg_data>`*`)`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_data.md)
  :

  Convert `eeg_data` to `data.frame`

- [`as.data.frame(`*`<eeg_epochs>`*`)`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_epochs.md)
  :

  Convert `eeg_epochs` object to data.frame

- [`as.data.frame(`*`<eeg_evoked>`*`)`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_evoked.md)
  :

  Convert `eeg_evoked` object to data frame

- [`as.data.frame(`*`<eeg_lm>`*`)`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_lm.md)
  :

  Convert `eeg_lm` to `data.frame`

- [`as.data.frame(`*`<eeg_stats>`*`)`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_stats.md)
  :

  Convert `eeg_stats` objects to data frames

- [`as.data.frame(`*`<eeg_tfr>`*`)`](https://craddm.github.io/eegUtils/reference/as.data.frame.eeg_tfr.md)
  :

  Convert `eeg_tfr` objects to a `data.frame`

## Accessors

Functions for accessing and modifying specific elements of an object

- [`channels()`](https://craddm.github.io/eegUtils/reference/channels.md)
  [`` `channels<-`() ``](https://craddm.github.io/eegUtils/reference/channels.md)
  : Modify channel information

- [`epochs()`](https://craddm.github.io/eegUtils/reference/epochs.md)
  [`` `epochs<-`() ``](https://craddm.github.io/eegUtils/reference/epochs.md)
  : Modify the epochs structure

- [`events()`](https://craddm.github.io/eegUtils/reference/events.md)
  [`` `events<-`() ``](https://craddm.github.io/eegUtils/reference/events.md)
  : Modify events structure

- [`channel_names()`](https://craddm.github.io/eegUtils/reference/channel_names.md)
  : Retrieve signal/channel names

- [`get_participant_id()`](https://craddm.github.io/eegUtils/reference/get_participant_id.md)
  [`get_recording()`](https://craddm.github.io/eegUtils/reference/get_participant_id.md)
  [`set_participant_id()`](https://craddm.github.io/eegUtils/reference/get_participant_id.md)
  [`set_recording()`](https://craddm.github.io/eegUtils/reference/get_participant_id.md)
  :

  Query and set elements of the `epochs` metadata structures

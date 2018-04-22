# eegUtils 0.1.15.dev

<a href="http://www.repostatus.org/#wip"><img src="http://www.repostatus.org/badges/latest/wip.svg" alt="Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public." /></a> [![Coverage Status](https://img.shields.io/codecov/c/github/craddm/eegUtils/master.svg)](https://codecov.io/github/craddm/eegUtils?branch=master) [![Build Status](https://travis-ci.org/craddm/eegUtils.svg?branch=master)](https://travis-ci.org/craddm/eegUtils)

Some helper utilities for plotting and processing EEG data in R. The package is in the early stages of development, and may be subject to a lot of changes.

Use devtools::install_github('craddm/eegUtils') to install it.

## Pre-processing functions
* import_raw() - for import BDF/EDF (BioSemi/European Data Format) EEG files 
* iir_filter() - for performing IIR Butterworth filtering on data  (FIR filtering in development)
* reref_eeg() - for re-referencing data
* epoch_data() - for creating epochs around trigger events
* tag_events() - for labelling events 
* interp_elecs() - spherical spline interpolation of EEG channels
* eeg_downsample() - filter and downsample data to a lower sampling rate.

## Plotting functions 
* topoplot() - plotting of topographies 
* plot_timecourse()/plot_butterfly()- plotting individual timecourses from electrodes or plotting all electrodes at once
* erp_scalp() - plotting ERP plots for individual electrodes in a topographical layout - thanks to Matti Vuorre!
* interactive_scalp() - a Shiny version of erp_scalp() that allows you to zoom in on specific electrodes.
* browse_data() - a Shiny gadget for interactively scrolling through EEG data (continous or epoched)

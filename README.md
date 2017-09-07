# eegUtils
<a href="http://www.repostatus.org/#wip"><img src="http://www.repostatus.org/badges/latest/wip.svg" alt="Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public." /></a>

Some helper utilities for plotting EEG data in R. The package is in the early stages of development, and may be subject to a lot of changes.

devtools::install_github('craddm/eegUtils') should install it.

There are currently a few working commands - 

* topoplot()
* plot_timecourse()/plot_butterfly() 
* erp_scalp() - thanks to Matti Vuorre!
* interactive_scalp() - a Shiny version of erp_scalp() that allows you to zoom in on specific electrodes.

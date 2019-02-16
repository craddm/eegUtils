.onAttach <- function(...) {

  packageStartupMessage(paste("Make sure to check for the latest development version at https://github.com/craddm/eegUtils!"))

}

utils::globalVariables(c("time",
                         "amplitude",
                         "electrode",
                         "epoch",
                         "cond_label",
                         "frequency",
                         "power",
                         "x",
                         "y",
                         "smooth_time",
                         ".",
                         "smooth_amp",
                         "pt",
                         "bl",
                         "..level..",
                         "condition",
                         "conditions",
                         "topos",
                         ".electrode",
                         "phi",
                         "theta",
                         "event_onset",
                         "variance",
                         "means",
                         "add_markers",
                         "plotlyOutput",
                         "value",
                         "var",
                         "measure",
                         "keep_rows",
                         "xyz_coords",
                         "component",
                         "signals"))

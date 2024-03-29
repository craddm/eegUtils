% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signal_decomposition.R
\name{eeg_decompose}
\alias{eeg_decompose}
\alias{eeg_decompose.eeg_epochs}
\title{Generalized eigenvalue decomposition based methods for EEG data}
\usage{
eeg_decompose(data, ...)

\method{eeg_decompose}{eeg_epochs}(
  data,
  sig_range,
  noise_range,
  method = "ssd",
  verbose = TRUE,
  order = 2,
  ...
)
}
\arguments{
\item{data}{An \code{eeg_data} object}

\item{...}{Additional parameters}

\item{sig_range}{Vector with two inputs, the lower and upper bounds of the frequency range of interest}

\item{noise_range}{Range of frequencies to be considered noise (e.g. bounds of flanker frequencies)}

\item{method}{Type of decomposition to apply. Currently only "ssd" is supported.}

\item{verbose}{Informative messages printed to console. Defaults to TRUE.}

\item{order}{Filter order for filter applied to signal/noise}
}
\value{
An \code{eeg_ICA} object. Note that
}
\description{
Implements a selection of Generalized Eigenvalue based decomposition methods
for EEG signals. Intended for isolating oscillations at specified
frequencies, decomposing channel-based data into components
reflecting distinct or combinations of sources of oscillatory signals.
Currently, spatio-spectral decomposition (Nikulin, Nolte, & Curio, 2011) and
Rhythmic Entrainment Source Separation (Cohen & Gulbinate, 2017) are
implemented. The key difference between the two is that the former returns
the results of the data-derived spatial filters applied to the
bandpass-filtered "signal" data, whereas the latter returns the results of the
filters applied to the original, broadband data.
}
\section{Methods (by class)}{
\itemize{
\item \code{eeg_epochs}: method for \code{eeg_epochs} objects
}}

\examples{
# The default method is Spatio-Spectral Decomposition, which returns
# spatially and temporally filtered source timecourses.
 decomposed <-
   eeg_decompose(demo_epochs,
                 sig_range = c(9, 11),
                 noise_range = c(8, 12),
                 method = "ssd")
 plot_psd(decomposed)
 # We can plot the spatial filters using `topoplot()`
 topoplot(decomposed, 1:2)
 plot_timecourse(decomposed, 1)
# method = "ress" returns spatially but not temporally filtered timecourses.
 with_RESS <-
   eeg_decompose(demo_epochs,
                 sig_range = c(9, 11),
                 noise_range = c(8, 12),
                 method = "ress")
 plot_psd(with_RESS)
 # The topographical plots are identical to those using "ssd", as the
 # spatial filters are the same.
 topoplot(with_RESS, 1:2)
 plot_timecourse(with_RESS, 1)
}
\references{
Cohen, M. X., & Gulbinate, R. (2017). Rhythmic entrainment source
separation: Optimizing analyses of neural responses to rhythmic sensory
stimulation. NeuroImage, 147, 43-56.
https://doi.org/10.1016/j.neuroimage.2016.11.036

Haufe, S., Dähne, S., & Nikulin, V. V. (2014). Dimensionality reduction for
the analysis of brain oscillations. NeuroImage, 101, 583–597.
https://doi.org/10.1016/j.neuroimage.2014.06.073

Nikulin, V. V., Nolte, G., & Curio, G. (2011). A novel method for reliable
and fast extraction of neuronal EEG/MEG oscillations on the basis of
spatio-spectral decomposition. NeuroImage, 55(4), 1528–1535.
https://doi.org/10.1016/j.neuroimage.2011.01.057
}
\seealso{
Other decompositions: 
\code{\link{run_ICA}()}
}
\author{
Matt Craddock \email{matt@mattcraddock.com}
}
\concept{decompositions}

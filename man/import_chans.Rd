% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/channel_management.R
\name{import_chans}
\alias{import_chans}
\title{Import channel locations from various file formats}
\usage{
import_chans(file_name, format = "spherical", file_format = "auto")
}
\arguments{
\item{file_name}{Name and full path of file to be loaded.}

\item{format}{If the file is not \code{.elc} format, "spherical", "geographic".
Default is "spherical".}

\item{file_format}{Default is \code{auto}, which will use the file extension to
determine file format. Other options include \code{ced}, \code{besa}, \code{elp}, \code{elc}}
}
\value{
A \code{tibble} containing electrode names and locations in several
different coordinate systems.
}
\description{
Currently only ASA \code{.elc} format with Cartesian x-y-z coordinates is
supported.
}
\author{
Matt Craddock \email{matt@mattcraddock.com}
}

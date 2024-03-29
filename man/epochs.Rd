% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epoch_queries.R
\name{epochs}
\alias{epochs}
\alias{epochs<-}
\title{Modify the epochs structure}
\usage{
epochs(data)

epochs(data) <- value
}
\arguments{
\item{data}{\code{eegUtils} object to view}

\item{value}{Structure to replace \code{epochs} structure with.}
}
\description{
Get or set the epochs structure of an \code{eegUtils} object.
}
\examples{
  epochs(demo_spatial)
}
\author{
Matt Craddock \email{matt@mattcraddock.com}
}

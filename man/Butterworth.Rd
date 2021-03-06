% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{Butterworth}
\alias{Butterworth}
\alias{Butterworth.numeric}
\alias{Butterworth.data.frame}
\alias{Butterworth.list}
\title{Butterworth Filter}
\usage{
Butterworth(x, ...)

\method{Butterworth}{numeric}(x, samplingRate, order, cutoffs, type, ...)

\method{Butterworth}{data.frame}(x, samplingRate, order, cutoffs, type, ...)

\method{Butterworth}{list}(x, samplingRate, order, cutoffs, type, ...)
}
\arguments{
\item{x}{A numeric vector or a data frame of numeric vectors.}

\item{...}{Placeholder for future methods.}

\item{samplingRate}{An atomic numeric, the sampling rate of \code{x}, in
hertz.}

\item{order}{An atomic numeric, the order of the filter.}

\item{cutoffs}{A numeric vector, the critical frequencies of the filter, in
hertz. For a low-pass or a high-pass filter, \code{cutoffs} should be of
length 1, denoting the cutoff frequency of the filter. For a pass-band or
stop-band filter, \code{cutoffs} should be of length 2, denoting the
corner frequencies of the band.}

\item{type}{A character string, either \code{"low"}, \code{"high"},
\code{"pass"}, or \code{"stop"}.}
}
\value{
The filtered data.
}
\description{
Apply a Butterworth filter to data.
}
\details{
\code{Butterworth} wraps \code{signal::butter} and \code{signal::filtfilt},
putting data arguments first so that filtering plays nice with the forward
pipe.  The filtering is done in the forward and reverse direction in order
to remove phase distortion from one-pass filtering.
}
\seealso{
\code{\link[signal]{butter}}, \code{\link[signal]{filtfilt}}
}

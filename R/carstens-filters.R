

#' Post-processing Filters
#'
#' Pre-defined finite impulse response (FIR) low-pass filters that are packaged
#' with the Carstens software on the control server. These filters are applied
#' to the recorded data by the Calcpos and Normpos applications during
#' post-processing. Each filter is represented as a Kaiser window.
#'
#' @format A data frame with 142 rows and 7 variables:
#' \itemize{
#'   \item \code{Name}: the name of the filter, as referenced by Calcpos or
#'         Normpos
#'   \item \code{PassbandLimit}: the (upper) edge of the filter's passband,
#'         in hertz
#'   \item \code{StopbandLimit}: the (lower) edge of the filter's stopband,
#'         in hertz
#'   \item \code{Attenuation}: the attenuation achieved within the stopband,
#'         in decibels
#'   \item \code{SamplingRate}: the sampling rate of the window function that
#'         implements the filter
#'   \item \code{DifferenceOrder}: for a given passband limit, stopband limit,
#'         attenuation, and sampling rate, the difference order of the window
#'         function that implements the filter: 0 denotes "position", 1 denotes
#'         "velocity", and 2 denotes "acceleration"
#'   \item \code{Window}: the values of the window function that implements the
#'         filter
#' }
#'
#' @source \code{/opt/age/bin/filter/} on Carstens control server.
#'
#' @export
CarstensFilters <- function() {
  tibble::as_tibble(carstens_filters)
}

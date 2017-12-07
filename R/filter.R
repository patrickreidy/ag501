

# Butterworth          #########################################################
################################################################################

#' S4 Generic: Butterworth
#'
#' Apply a Butterworth filter to data.
#'
#' \code{Butterworth} wraps \code{signal::butter} and \code{signal::filtfilt},
#' putting data arguments first so that filtering plays nice with the forward
#' pipe.  The filtering is done in the forward and reverse direction in order
#' to remove phase distortion from one-pass filtering.
#'
#' @param x A numeric vector or a data frame of numeric vectors.
#'
#' @return The filtered data.
#'
#' @seealso \code{\link[signal]{butter}}, \code{\link[signal]{filtfilt}}
#'
#' @name Butterworth
#' @export
methods::setGeneric(
  name = "Butterworth",
  def = function(x, ...) {
    standardGeneric("Butterworth")
  }
)


#' @param samplingRate An atomic numeric, the sampling rate of \code{x}, in
#'   hertz.
#' @param order An atomic numeric, the order of the filter.
#' @param cutoffs A numeric vector, the critical frequencies of the filter, in
#'   hertz. For a low-pass or a high-pass filter, \code{cutoffs} should be of
#'   length 1, denoting the cutoff frequency of the filter. For a pass-band or
#'   stop-band filter, \code{cutoffs} should be of length 2, denoting the
#'   corner frequencies of the band.
#' @param type A character string, either \code{"low"}, \code{"high"},
#'   \code{"pass"}, or \code{"stop"}. Default is \code{"low"}.
#'
#' @usage \S4method{Butterworth}{numeric}(x, samplingRate, order, cutoffs, type)
#'
#' @name Butterworth,numeric-method
#' @rdname Butterworth
methods::setMethod(
  f = "Butterworth",
  signature = c(x = "numeric"),
  definition = function(x, samplingRate, order, cutoffs, type = "low") {
    signal::butter(n = order, W = cutoffs/samplingRate/2, type = type) %>%
      signal::filtfilt(x)
  }
)


#' @param suffix A character string. If \code{x} is a data frame, the names of the
#'   returned data frame are equal to the names of \code{x}, each suffixed by
#'   \code{suffix}. Default is \code{""}.
#'
#' @usage \S4method{Butterworth}{data.frame}(x, samplingRate, order, cutoffs, type, suffix)
#'
#' @name Butterworth,data.frame-method
#' @rdname Butterworth
methods::setMethod(
  f = "Butterworth",
  signature = c(x = "data.frame"),
  definition = function(x, samplingRate, order, cutoffs, type = "low", suffix = "") {
    .filtered <- purrr::map_df(x, Butterworth,
                               samplingRate = samplingRate, order = order,
                               cutoffs = cutoffs, type = type)
    .named <- purrr::set_names(.filtered, stringr::str_c(names(.filtered), suffix))
    return(.named)
  }
)




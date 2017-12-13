

# Butterworth          #########################################################
################################################################################

#' Butterworth Filter
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
#' @export
Butterworth <- function(x, ...) {
  UseMethod("Butterworth", x)
}


#' @param samplingRate An atomic numeric, the sampling rate of \code{x}, in
#'   hertz.
#' @param order An atomic numeric, the order of the filter.
#' @param cutoffs A numeric vector, the critical frequencies of the filter, in
#'   hertz. For a low-pass or a high-pass filter, \code{cutoffs} should be of
#'   length 1, denoting the cutoff frequency of the filter. For a pass-band or
#'   stop-band filter, \code{cutoffs} should be of length 2, denoting the
#'   corner frequencies of the band.
#' @param type A character string, either \code{"low"}, \code{"high"},
#'   \code{"pass"}, or \code{"stop"}.
#' @param ... Placeholder for future methods.
#'
#'
#' @rdname Butterworth
#' @export
Butterworth.numeric <- function(x, samplingRate, order, cutoffs, type, ...) {
  signal::butter(n = order, W = cutoffs/samplingRate/2, type = type) %>%
    signal::filtfilt(x)
}


#' @rdname Butterworth
#' @export
Butterworth.data.frame <- function(x, samplingRate, order, cutoffs, type, ...) {
  .channel_data <- dplyr::select(x, -.data$Time)
  .filtered_data <- purrr::map_df(.channel_data, Butterworth,
                                  samplingRate = samplingRate, order = order,
                                  cutoffs = cutoffs, type = type)
  .position_data <- dplyr::bind_cols(dplyr::select(x, .data$Time), .filtered_data)
  return(.position_data)
}


#' @rdname Butterworth
#' @export
Butterworth.list <- function(x, samplingRate, order, cutoffs, type, ...) {
  if (length(samplingRate) != length(x) & length(samplingRate) != 1) {
    .stop_msg <- sprintf("length of @samplingRate (%d) must be either length of @x (%d) or 1.",
                         length(samplingRate), length(x))
    stop(.stop_msg)
  } else if (length(samplingRate) == 1) {
    samplingRate <- rep(samplingRate, length(x))
  }
  if (length(order) != length(x) & length(order) != 1) {
    .stop_msg <- sprintf("length of @order (%d) must be either length of @x (%d) or 1.",
                         length(order), length(x))
    stop(.stop_msg)
  } else if (length(order) == 1) {
    order <- rep(order, length(x))
  }
  if (length(cutoffs) != length(x) & length(cutoffs) != 1) {
    .stop_msg <- sprintf("length of @cutoffs (%d) must be either length of @x (%d) or 1.",
                         length(cutoffs), length(x))
    stop(.stop_msg)
  } else if (length(cutoffs) == 1) {
    cutoffs <- rep(cutoffs, length(x))
  }
  if (length(type) != length(x) & length(type) != 1) {
    .stop_msg <- sprintf("length of @type (%d) must be either length of @x (%d) or 1.",
                         length(type), length(x))
    stop(.stop_msg)
  } else if (length(type) == 1) {
    type <- rep(type, length(x))
  }
  .filtered <-
    purrr::pmap(list(x, samplingRate, order, cutoffs, type),
                function(.d, .sr, .or, .ct, .tp) {
                  Butterworth(x = .d, samplingRate = .sr, order = .or,
                              cutoffs = .ct, type = .tp)
                })
  return(.filtered)
}


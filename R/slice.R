

#' Time Slice
#'
#' Slice a contiguous temporal interval of kinematic data.
#'
#' @param x A data.frame or a list of data.frames that contain kinematic data,
#'   structured such that there is a column named \code{Time} and other columns
#'   that record the channels of Carstens sensors.
#' @param from A numeric, the time from which to slice.
#' @param to A numeric, the time to which to slice.
#' @param ... Placeholder for future methods.
#'
#' @return An object that has the same class as \code{x}. If \code{x} is a
#'   data.frame, then the returned data.frame is like \code{x}, but filtered
#'   to just those rows whose value of \code{Time} falls within the
#'   \code{[from, to]} interval. If \code{x} is a list, then a list of
#'   time-sliced data.frames is returned.
#'
#' @export
TimeSlice <- function(x, from, to, ...) {
  UseMethod("TimeSlice", x)
}


#' @rdname TimeSlice
#' @export
TimeSlice.data.frame <- function(x, from, to, ...) {
  .sliced <- dplyr::filter(x, from <= .data$Time, .data$Time <= to)
  return(.sliced)
}


#' @rdname TimeSlice
#' @export
TimeSlice.list <- function(x, from, to, ...) {
  if (length(from) == 1) {
    from <- rep(from, times = length(x))
  } else if (length(from) != length(x)) {
    .stop_msg <- sprintf("length of @from (%d) must be either length of @x (%d) or 1.",
                         length(from), length(x))
    stop(.stop_msg)
  }
  if (length(to) == 1) {
    to <- rep(to, times = length(x))
  } else if (length(to) != length(x)) {
    .stop_msg <- sprintf("length of @to (%d) must be either length of @x (%d) or 1.",
                         length(to), length(x))
    stop(.stop_msg)
  }
  .sliced <- purrr::pmap(list(x, from, to), function(.x, .f, .t) {
    TimeSlice(x = .x, from = .f, to = .t)
  })
  return(.sliced)
}

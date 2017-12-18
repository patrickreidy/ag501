

# Find Maxima          #########################################################
################################################################################

#' Find Local Maxima
#'
#' Find the local maxima within one or more channels of kinematic data.
#'
#' @param x A data.frame (or a list of data.frames) that contain(s) kinematic
#'   data (e.g., position, velocity, speed, or acceleration data). \code{x}
#'   must have one column named \code{Time} and at least one other column.
#' @param channels A character vector, each element of which should match
#'   exactly the name of some column in \code{x} other than \code{Time}.
#'   Default matches all columns of \code{x} other than \code{Time}.
#' @param matches A pattern used to match columns in \code{x}. A non-null
#'   value for this argument will supercede a non-null value for the
#'   \code{channels} argument. Default matches all columns of \code{x} other
#'   than \code{Time}.
#' @param from,to Atomic numerics. Local maxima will only be found within the
#'   \code{[from,to]} interval. Default values will include all available
#'   data in the channel.
#' @param alpha A proportion within the range \code{[0,1)}. Local maxima in a
#'   given channel will only be found if their value exceeds \code{100*alpha}
#'   percent of the data range in that channel.
#' @param ... Placeholder for future methods.
#'
#' @return A data table with the following variables:
#'   \itemize{
#'     \item \code{Channel}: the name of the channel in which the local
#'       maximum occurs.
#'     \item \code{Time}: the time at which the local occurs.
#'     \item \code{Value}: the amplitude value of the local maximum.
#'   }
#'
#' @export
FindMaxima <- function(x, ...) {
  UseMethod("FindMaxima", x)
}


#' @rdname FindMaxima
#' @export
FindMaxima.data.frame <- function(x, channels = NULL, matches = NULL, from = -Inf, to = Inf, alpha = 0.1, ...) {
  .FindMaximaInSingleChannel <- function(.channel_data, .a = alpha) {
    .channel_name <- setdiff(names(.channel_data), "Time")
    .channel_values <- dplyr::pull(.channel_data, rlang::UQ(.channel_name))
    .value_thresh <- min(.channel_values, na.rm = TRUE) +
      (.a * (max(.channel_values, na.rm = TRUE) - min(.channel_values, na.rm = TRUE)))
    .local_maxs <-
      purrr::set_names(.channel_data, nm = c("Time", "Value")) %>%
      dplyr::mutate(Channel = .channel_name) %>%
      dplyr::mutate(IncreasingBefore = (.data$Value - dplyr::lag(.data$Value)) > 0) %>%
      dplyr::mutate(DecreasingAfter = (.data$Value - dplyr::lead(.data$Value)) > 0) %>%
      dplyr::filter(.data$IncreasingBefore, .data$DecreasingAfter, .data$Value >= .value_thresh) %>%
      dplyr::select(.data$Channel, .data$Time, .data$Value)
    return(.local_maxs)
  }
  if (!is.null(matches)) {
    .channels <- grep(pattern = matches, x = setdiff(names(x), "Time"), value = TRUE)
  } else if (!is.null(channels)) {
    .channels <- setdiff(channels, "Time")
  } else {
    .channels <- setdiff(names(x), "Time")
  }
  .channels_data <-
    purrr::map(.channels, function(.ch) {
      dplyr::select(x, .data$Time, rlang::UQ(.ch)) %>%
        dplyr::filter(from <= .data$Time, .data$Time <= to)
    })
  .local_maxs <-
    purrr::map(.channels_data, .FindMaximaInSingleChannel, .a = alpha)
  .local_maxs <- purrr::reduce(.local_maxs, dplyr::bind_rows)
  return(.local_maxs)
}


#' @rdname FindMaxima
#' @export
FindMaxima.list <- function(x, channels = NULL, matches = NULL, from = -Inf, to = Inf, alpha = 0.1, ...) {
  if (!is.list(matches)) {
    matches <- list(matches)
  }
  if (length(matches) == 1) {
    matches <- rep(matches, times = length(x))
  } else if (length(matches) != length(x)) {
    .stop_msg <- sprintf("length of @matches (%d) must be either length of @x (%d) or 1.",
                         length(matches), length(x))
    stop(.stop_msg)
  }
  if (!is.list(channels)) {
    channels <- list(channels)
  }
  if (length(channels) == 1) {
    channels <- rep(channels, times = length(x))
  } else if (length(matches) != length(x)) {
    .stop_msg <- sprintf("length of @channels (%d) must be either length of @x (%d) or 1.",
                         length(channels), length(x))
    stop(.stop_msg)
  }
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
  if (length(alpha) == 1) {
    alpha = rep(alpha, times = length(x))
  } else if (length(alpha) != length(x)) {
    .stop_msg <- sprintf("length of @alpha (%d) must be either length of @x (%d) or 1.",
                         length(alpha), length(x))
    stop(.stop_msg)
  }
  .local_maxs <- purrr::pmap(list(x, channels, matches, from, to, alpha),
                             function(.x, .c, .m, .f, .t, .a) {
                               FindMaxima(.x, channels = .c, matches = .m,
                                          from = .f, to = .t, alpha = .a)
                             })
  return(.local_maxs)
}










# Find Minima          #########################################################
################################################################################

#' Find Local Minima
#'
#' Find the local minima within one or more channels of kinematic data.
#'
#' @param x A data.frame (or a list of data.frames) that contain(s) kinematic
#'   data (e.g., position, velocity, speed, or acceleration data). \code{x}
#'   must have one column named \code{Time} and at least one other column.
#' @param channels A character vector, each element of which should match
#'   exactly the name of some column in \code{x} other than \code{Time}.
#'   Default matches all columns of \code{x} other than \code{Time}.
#' @param matches A pattern used to match columns in \code{x}. A non-null
#'   value for this argument will supercede a non-null value for the
#'   \code{channels} argument. Default matches all columns of \code{x} other
#'   than \code{Time}.
#' @param from,to Atomic numerics. Local maxima will only be found within the
#'   \code{[from,to]} interval. Default values will include all available
#'   data in the channel.
#' @param ... Placeholder for future methods.
#'
#' @return A data table with the following variables:
#'   \itemize{
#'     \item \code{Channel}: the name of the channel in which the local
#'       maximum occurs.
#'     \item \code{Time}: the time at which the local occurs.
#'     \item \code{Value}: the amplitude value of the local minimum.
#'   }
#'
#' @export
FindMinima <- function(x, ...) {
  UseMethod("FindMinima", x)
}


#' @rdname FindMinima
#' @export
FindMinima.data.frame <- function(x, channels = NULL, matches = NULL, from = -Inf, to = Inf, ...) {
  .FindMinimaInSingleChannel <- function(.channel_data) {
    .channel_name <- setdiff(names(.channel_data), "Time")
    .local_mins <-
      purrr::set_names(.channel_data, nm = c("Time", "Value")) %>%
      dplyr::mutate(Channel = .channel_name) %>%
      dplyr::mutate(DecreasingBefore = (.data$Value - dplyr::lag(.data$Value)) < 0) %>%
      dplyr::mutate(IncreasingAfter = (.data$Value - dplyr::lead(.data$Value)) < 0) %>%
      dplyr::filter(.data$DecreasingBefore, .data$IncreasingAfter) %>%
      dplyr::select(.data$Channel, .data$Time, .data$Value)
    return(.local_mins)
  }
  if (!is.null(matches)) {
    .channels <- grep(pattern = matches, x = setdiff(names(x), "Time"), value = TRUE)
  } else if (!is.null(channels)) {
    .channels <- setdiff(channels, "Time")
  } else {
    .channels <- setdiff(names(x), "Time")
  }
  .channels_data <-
    purrr::map(.channels, function(.ch) {
      dplyr::select(x, .data$Time, rlang::UQ(.ch)) %>%
        dplyr::filter(from <= .data$Time, .data$Time <= to)
    })
  .local_mins <-
    purrr::map(.channels_data, .FindMinimaInSingleChannel)
  .local_mins <- purrr::reduce(.local_mins, dplyr::bind_rows)
  return(.local_mins)
}


#' @rdname FindMinima
#' @export
FindMinima.list <- function(x, channels = NULL, matches = NULL, from = -Inf, to = Inf, ...) {
  if (!is.list(matches)) {
    matches <- list(matches)
  }
  if (length(matches) == 1) {
    matches <- rep(matches, times = length(x))
  } else if (length(matches) != length(x)) {
    .stop_msg <- sprintf("length of @matches (%d) must be either length of @x (%d) or 1.",
                         length(matches), length(x))
    stop(.stop_msg)
  }
  if (!is.list(channels)) {
    channels <- list(channels)
  }
  if (length(channels) == 1) {
    channels <- rep(channels, times = length(x))
  } else if (length(matches) != length(x)) {
    .stop_msg <- sprintf("length of @channels (%d) must be either length of @x (%d) or 1.",
                         length(channels), length(x))
    stop(.stop_msg)
  }
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
  .local_mins <- purrr::pmap(list(x, channels, matches, from, to),
                             function(.x, .c, .m, .f, .t) {
                               FindMinima(.x, channels = .c, matches = .m,
                                          from = .f, to = .t)
                             })
  return(.local_mins)
}










# Find Landmarks          ######################################################
################################################################################

#' Find Gestural Landmarks
#'
#' Find the gestural onset, target onset, target offset, and gestural offset
#' within a speed time-series.
#'
#' Gestural landmarks may be found only within a time series that denotes the
#' speed of a sensor (or of a derived quantity, such as lip aperture) in one
#' or more dimensions. Within a speed time-series, a gesture is comprises three
#' local minima (\code{m1}, \code{m3}, \code{m5}), with a local maximum
#' (\code{m2}, \code{m4}) between each local minima: i.e., these local minima
#' and maxima are ordered such that \code{m1 < m2 < m3 < m4 < m5}. The first
#' half of the gesture, during which the constriction is formed, corresponds to
#' the interval \code{[m1, m3)}, where the speed values trace a concave downward
#' trajectory. The second half of the gesture, during which the constriction
#' is released, corresponds to the interval \code{[m3, m5]}, where again the
#' speed values trace a concave downward trajectory.
#'
#' A gestural landmark corresponds to the point in time between two consecutive
#' local extrema, where the speed time-series is equal to some proportion
#' (\code{theta}) of the range between the two local extrema. For example, if
#' \code{theta = 0.2}, then the gestural onset is the point in time where the
#' speed is equal to 20 percent of the range between the speed at \code{m1},
#' and the speed at \code{m2}: \code{speed(m1) + (0.2 * (speed(m2) - speed(m1)))}.
#'
#' In order to find landmarks, the user must provide a fixed point in time that
#' is near to the local minima of the onset of the gesture (i.e., via
#' the \code{onsetNear} argument) or the target (i.e., via the \code{targetNear}
#' argument).
#'
#' @param x A data.frame (or a list of data.frames) that contain(s) kinematic
#'   data (e.g., position, velocity, speed, or acceleration data). \code{x}
#'   must have one column named \code{Time} and at least one other column.
#' @param channels A character vector, each element of which should match
#'   exactly the name of some column in \code{x} other than \code{Time}.
#'   Default matches all columns of \code{x} other than \code{Time}.
#' @param matches A pattern used to match columns in \code{x}. A non-null
#'   value for this argument will supercede a non-null value for the
#'   \code{channels} argument. Default matches all columns of \code{x} other
#'   than \code{Time}.
#' @param from,to Atomic numerics. Local maxima will only be found within the
#'   \code{[from,to]} interval. Default values will include all available
#'   data in the channel.
#' @param onsetNear,targetNear Atomic numerics.
#' @param theta A proportion within the range \code{[0,1)}. Given consecutive
#'   local maximum \code{max} and local minimum \code{min}, the gestural
#'   landmark \code{mark} is the point in time between the local maximum
#'   and minimum such that:
#'   \code{value(mark) = value(min) + (theta * (value(max) - value(min)))}.
#'   Default is \code{0.2}, which is the "industry standard".
#' @param alpha A proportion within the range \code{[0,1)}. Local maxima in a
#'   given channel will only be found if their value exceeds \code{100*alpha}
#'   percent of the data range in that channel.
#' @param ... Placeholder for future methods.
#'
#' @return A data table with the following variables:
#'   \itemize{
#'     \item \code{Channel}: the name of the channel in which the local
#'       maximum occurs.
#'     \item \code{Landmark}: the name of the landmark: \code{GestureOnset},
#'       \code{TargetOnset}, \code{TargetOffset}, \code{GestureOffset}.
#'     \item \code{Time}: the time at which the local occurs.
#'     \item \code{Value}: the amplitude value of the local minimum.
#'   }
#'
#' @export
FindLandmarks <- function(x, ...) {
  UseMethod("FindLandmarks", x)
}

#' @rdname FindLandmarks
#' @export
FindLandmarks.data.frame <- function(x, channels = NULL, matches = NULL, from = -Inf, to = Inf, onsetNear = NULL, targetNear = NULL, theta = 0.2, alpha = 0.1, ...) {
  .FindLandmarksInSingleChannel <- function(.channel_data) {
    .maxima <- FindMaxima(x = .channel_data, from = from, to = to, alpha = alpha)
    .minima <- FindMinima(x = .channel_data, from = from, to = to)
    if (!is.null(targetNear)) {
      .i <- which.min(abs(.minima$Time - targetNear))
      .m1_i <- .i - 1
    } else if (!is.null(onsetNear)) {
      .m1_i <- which.min(abs(.minima$Time - onsetNear))
    } else {
      stop("either @onsetNear or @targetNear must be a non-NULL numeric")
    }
    .m1 <- .minima[.m1_i, ]
    .m2 <- .maxima %>% dplyr::filter(.data$Time > .m1$Time) %>% dplyr::slice(1)
    .m3 <- .minima[.m1_i+1, ]
    .m4 <- .maxima %>% dplyr::filter(.data$Time > .m3$Time) %>% dplyr::slice(1)
    .m5 <- .minima[.m1_i+2, ]
    .Landmark <- function(.min, .max, .landmark) {
      .channel_name <- setdiff(names(.channel_data), "Time")
      .earlier_extrema <- min(c(.min$Time, .max$Time))
      .later_extrema <- max(c(.min$Time, .max$Time))
      .threshold_value <- .min$Value + (theta * (.max$Value - .min$Value))
      .landmark_data <-
        .channel_data %>%
        purrr::set_names(nm = c("Time", "Value")) %>%
        dplyr::filter(.earlier_extrema <= .data$Time, .data$Time <= .later_extrema) %>%
        dplyr::arrange(abs(.data$Value - .threshold_value)) %>%
        dplyr::slice(1) %>%
        dplyr::mutate(Channel = .channel_name, Landmark = .landmark) %>%
        dplyr::select(.data$Channel, .data$Landmark, .data$Time, .data$Value)
      return(.landmark_data)
    }
    .mins <- list(.m1, .m3, .m3, .m5)
    .maxs <- list(.m2, .m2, .m4, .m4)
    .marks <- list("GestureOnset", "TargetOnset", "TargetOffset", "GestureOffset")
    .landmarks <-
      purrr::pmap(list(.mins, .maxs, .marks),
                  function(.n, .x, .k) {.Landmark(.min = .n, .max = .x, .landmark = .k)}) %>%
      purrr::reduce(dplyr::bind_rows)
    return(.landmarks)
  }
  if (!is.null(matches)) {
    .channels <- grep(pattern = matches, x = setdiff(names(x), "Time"), value = TRUE)
  } else if (!is.null(channels)) {
    .channels <- setdiff(channels, "Time")
  } else {
    .channels <- setdiff(names(x), "Time")
  }
  .channels_data <-
    purrr::map(.channels, function(.ch) {
      dplyr::select(x, .data$Time, rlang::UQ(.ch)) %>%
        dplyr::filter(from <= .data$Time, .data$Time <= to)
    })
  .landmarks <-
    purrr::map(.channels_data, .FindLandmarksInSingleChannel)
  .landmarks <- purrr::reduce(.landmarks, dplyr::bind_rows)
  return(.landmarks)
}


#' @rdname FindLandmarks
#' @export
FindLandmarks.list <- function(x, channels = NULL, matches = NULL, from = -Inf, to = Inf, onsetNear = NULL, targetNear = NULL, theta = 0.2, alpha = 0.1, ...) {
  if (!is.list(matches)) {
    matches <- list(matches)
  }
  if (!is.null(matches)) {
    if (length(matches) == 1) {
      matches <- rep(matches, times = length(x))
    } else if (length(matches) != length(x)) {
      .stop_msg <- sprintf("length of @matches (%d) must be either length of @x (%d) or 1.",
                           length(matches), length(x))
      stop(.stop_msg)
    }
  }
  if (!is.list(channels)) {
    channels <- list(channels)
  }
  if (length(channels) == 1) {
    channels <- rep(channels, times = length(x))
  } else if (length(matches) != length(x)) {
    .stop_msg <- sprintf("length of @channels (%d) must be either length of @x (%d) or 1.",
                         length(channels), length(x))
    stop(.stop_msg)
  }
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
  if (is.null(onsetNear)) {
    onsetNear <- list(onsetNear)
  }
  if (length(onsetNear) == 1) {
    onsetNear <- rep(onsetNear, times = length(x))
  } else if (length(onsetNear) != length(x)) {
    .stop_msg <- sprintf("length of @onsetNear (%d) must be either length of @x (%d) or 1.",
                         length(onsetNear), length(x))
    stop(.stop_msg)
  }
  if (is.null(targetNear)) {
    targetNear <- list(targetNear)
  }
  if (length(targetNear) == 1) {
    targetNear <- rep(targetNear, times = length(x))
  } else if (length(targetNear) != length(x)) {
    .stop_msg <- sprintf("length of @targetNear (%d) must be either length of @x (%d) or 1.",
                         length(targetNear), length(x))
    stop(.stop_msg)
  }
  if (length(theta) == 1) {
    theta = rep(theta, times = length(x))
  } else if (length(theta) != length(x)) {
    .stop_msg <- sprintf("length of @theta (%d) must be either length of @x (%d) or 1.",
                         length(theta), length(x))
    stop(.stop_msg)
  }
  if (length(alpha) == 1) {
    alpha = rep(alpha, times = length(x))
  } else if (length(alpha) != length(x)) {
    .stop_msg <- sprintf("length of @alpha (%d) must be either length of @x (%d) or 1.",
                         length(alpha), length(x))
    stop(.stop_msg)
  }
  .landmarks <- purrr::pmap(list(x, channels, matches, from, to, onsetNear, targetNear, theta, alpha),
                            function(.x, .c, .m, .f, .t, .on, .tn, .th, .al) {
                              FindLandmarks(.x, channels = .c, matches = .m,
                                            from = .f, to = .t,
                                            onsetNear = .on, targetNear = .tn,
                                            theta = .th, alpha = .al)
                            })
  return(.landmarks)
}

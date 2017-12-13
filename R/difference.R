

# Forward Difference          ##################################################
################################################################################

#' Forward Difference
#'
#' Given a numeric series \code{x}, compute its \code{n}-offset forward
#' difference \code{y}, such that \code{y[i] = x[i+n] - x[i]}.
#'
#' @param x A numeric vector or a data frame of numeric vectors.
#'
#' @return The forward difference of \code{x} if \code{x} is numeric, or of
#'   each column of \code{x} if \code{x} is a data frame.
#'
#' @export
ForwardDifference <- function(x, ...) {
  UseMethod("ForwardDifference", x)
}


#' @param n An atomic numeric, the number of samples forward from the current
#'   sample used in the forward difference. Default is \code{1}.
#' @param order An atomic numeric, the number of times the forward-difference
#'   operator is applied to \code{x}. Default is \code{1}.
#' @param samplingRate An atomic numeric, the sampling rate of \code{x}, in
#'   hertz. The raw forward difference of \code{x} is divided by the sampling
#'   rate so that the forward difference can be interepreted as an estimate of
#'   the derivative of \code{x}. Default is \code{n}, which yields the raw
#'   forward difference.
#' @param ... Placeholder for future methods.
#'
#' @rdname ForwardDifference
#' @export
ForwardDifference.numeric <- function(x, n = 1, order = 1, samplingRate = n, ...) {
  .x <- x
  .o <- 0
  while (.o < order) {
    .x <- (dplyr::lead(x = .x, n = n) - .x) / (n / samplingRate)
    .o <- .o + 1
  }
  return(.x)
}


#' @rdname ForwardDifference
#' @export
ForwardDifference.data.frame <- function(x, n = 1, order = 1, ...) {
  .x <- x
  .o <- 0
  while(.o < order) {
    .data <- dplyr::select(.x, -.data$Time)
    .dx <- purrr::map_df(.data, function(.d) {
      ForwardDifference(.d, n = n, order = 1) / ForwardDifference(.x$Time, n = n, order = 1)
    })
    .x <- dplyr::bind_cols(dplyr::select(.x, .data$Time), .dx)
    .o <- .o + 1
  }
  return(.x)
}


#' @rdname ForwardDifference
#' @export
ForwardDifference.list <- function(x, n = 1, order = 1, ...) {
  if (length(n) != length(x) & length(n) != 1) {
    .stop_msg <- sprintf("length of @n (%d) must be either length of @x (%d) or 1.",
                         length(n), length(x))
    stop(.stop_msg)
  } else if (length(n) == 1) {
    n <- rep(n, length(x))
  }
  if (length(order) != length(x) & length(order) != 1) {
    .stop_msg <- sprintf("length of @order (%d) must be either length of @x (%d) or 1.",
                         length(order), length(x))
    stop(.stop_msg)
  } else if (length(order) == 1) {
    order <- rep(order, length(x))
  }
  .differenced <-
    purrr::pmap(list(x, n, order), function(.x, .n, .o) {
      ForwardDifference(x = .x, n = .n, order = .o)
    })
  return(.differenced)
}






# Central Difference          ##################################################
################################################################################

#' Central Difference
#'
#' Given a numeric series \code{x}, compute its \code{n}-offset central
#' difference \code{y}, such that \code{y[i] = x[i+n] - x[i-n]}.
#'
#' @param x A numeric vector or a data frame of numeric vectors.
#'
#' @return The central difference of \code{x} if \code{x} is numeric, or of
#'   each column of \code{x} if \code{x} is a data frame.
#'
#' @export
CentralDifference <- function(x, ...) {
  UseMethod("CentralDifference", x)
}


#' @param n An atomic numeric, the number of samples offset from the current
#'   sample used in the central difference. Default is \code{1}.
#' @param order An atomic numeric, the number of times the central-difference
#'   operator is applied to \code{x}. Default is \code{1}.
#' @param samplingRate An atomic numeric, the sampling rate of \code{x}, in
#'   hertz. The raw central difference of \code{x} is divided by the sampling
#'   rate so that the central difference can be interepreted as an estimate of
#'   the derivative of \code{x}. Default is \code{2*n}, which yields the raw
#'   central difference.
#' @param ... Placeholder for future methods.
#'
#' @rdname CentralDifference
#' @export
CentralDifference.numeric <- function(x, n = 1, order = 1, samplingRate = 2*n, ...) {
  .x <- x
  .o <- 0
  while (.o < order) {
    .x <- (dplyr::lead(x = .x, n = n) - dplyr::lag(x = .x, n = n)) / (2*n / samplingRate)
    .o <- .o + 1
  }
  return(.x)
}


#' @rdname CentralDifference
#' @export
CentralDifference.data.frame <- function(x, n = 1, order = 1, ...) {
  .x <- x
  .o <- 0
  while(.o < order) {
    .data <- dplyr::select(.x, -.data$Time)
    .dx <- purrr::map_df(.data, function(.d) {
      CentralDifference(.d, n = n, order = 1) / CentralDifference(.x$Time, n = n, order = 1)
    })
    .x <- dplyr::bind_cols(dplyr::select(.x, .data$Time), .dx)
    .o <- .o + 1
  }
  return(.x)
}


#' @rdname ForwardDifference
#' @export
CentralDifference.list <- function(x, n = 1, order = 1, ...) {
  if (length(n) != length(x) & length(n) != 1) {
    .stop_msg <- sprintf("length of @n (%d) must be either length of @x (%d) or 1.",
                         length(n), length(x))
    stop(.stop_msg)
  } else if (length(n) == 1) {
    n <- rep(n, length(x))
  }
  if (length(order) != length(x) & length(order) != 1) {
    .stop_msg <- sprintf("length of @order (%d) must be either length of @x (%d) or 1.",
                         length(order), length(x))
    stop(.stop_msg)
  } else if (length(order) == 1) {
    order <- rep(order, length(x))
  }
  .differenced <-
    purrr::pmap(list(x, n, order), function(.x, .n, .o) {
      CentralDifference(x = .x, n = .n, order = .o)
    })
  return(.differenced)
}



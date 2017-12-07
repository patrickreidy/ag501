

# Forward Difference          ##################################################
################################################################################

#' S4 Generic: ForwardDifference
#'
#' Given a numeric series \code{x}, compute its \code{n}-offset forward
#' difference \code{y}, such that \code{y[i] = x[i+n] - x[i]}.
#'
#' @param x A numeric vector or a data frame of numeric vectors.
#'
#' @return The forward difference of \code{x} if \code{x} is numeric, or of
#'   each column of \code{x} if \code{x} is a data frame.
#'
#' @name ForwardDifference
#' @export
methods::setGeneric(
  name = "ForwardDifference",
  def = function(x, ...) {
    standardGeneric("ForwardDifference")
  }
)


#' @param samplingRate An atomic numeric, the sampling rate of \code{x}, in
#'   hertz. The raw forward difference of \code{x} is divided by the sampling
#'   rate so that the forward difference can be interepreted as an estimate of
#'   the derivative of \code{x}. Default is \code{1}.
#' @param n An atomic numeric, the number of samples forward from the current
#'   sample used in the forward difference. Default is \code{1}.
#' @param order An atomic numeric, the number of times the forward-difference
#'   operator is applied to \code{x}. Default is \code{1}.
#'
#' @usage \S4method{ForwardDifference}{numeric}(x, samplingRate, n, order)
#'
#' @name ForwardDifference,numeric-method
#' @rdname ForwardDifference
methods::setMethod(
  f = "ForwardDifference",
  signature = c(x = "numeric"),
  definition = function(x, samplingRate = 1, n = 1, order = 1) {
    .x <- x
    .o <- 0
    while (.o < order) {
      .x <- (dplyr::lead(x = .x, n = n) - .x) / (n / samplingRate)
      .o <- .o + 1
    }
    return(.x)

  }
)


#' @param suffix A character string. If \code{x} is a data frame, the names of the
#'   returned data frame are equal to the names of \code{x}, each suffixed by
#'   \code{suffix}. Default is \code{""}.
#'
#' @usage \S4method{ForwardDifference}{data.frame}(x, samplingRate, n, order, suffix)
#'
#' @name ForwardDifference,data.frame-method
#' @rdname ForwardDifference
methods::setMethod(
  f = "ForwardDifference",
  signature = c(x = "data.frame"),
  definition = function(x, samplingRate = 1, n = 1, order = 1, suffix = "") {
    .diffed <- purrr::map_df(x, ForwardDifference,
                             samplingRate = samplingRate, n = n)
    .named <- purrr::set_names(.diffed, stringr::str_c(names(.diffed), suffix))
    return(.named)
  }
)





# Central Difference          ##################################################
################################################################################

#' S4 Generic: CentralDifference
#'
#' Given a numeric series \code{x}, compute its \code{n}-offset central
#' difference \code{y}, such that \code{y[i] = x[i+n] - x[i-n]}.
#'
#' @param x A numeric vector or a data frame of numeric vectors.
#' @return The central difference of \code{x} if \code{x} is numeric, or of
#'   each column of \code{x} if \code{x} is a data frame.
#'
#' @name CentralDifference
#' @export
methods::setGeneric(
  name = "CentralDifference",
  def = function(x, ...) {
    standardGeneric("CentralDifference")
  }
)


#' @param samplingRate An atomic numeric, the sampling rate of \code{x}, in
#'   hertz. The raw central difference of \code{x} is divided by the sampling
#'   rate so that the central difference can be interepreted as an estimate of
#'   the derivative of \code{x}. Default is \code{1}.
#' @param n An atomic numeric, the number of samples offset from the current
#'   sample used in the central difference. Default is \code{1}.
#' @param order An atomic numeric, the number of times the central-difference
#'   operator is applied to \code{x}. Default is \code{1}.
#'
#' @usage \S4method{CentralDifference}{numeric}(x, samplingRate, n, order)
#'
#' @name CentralDifference,numeric-method
#' @rdname CentralDifference
methods::setMethod(
  f = "CentralDifference",
  signature = c(x = "numeric"),
  definition = function(x, samplingRate = 1, n = 1, order = 1) {
    .x <- x
    .o <- 0
    while (.o < order) {
      .x <- (dplyr::lead(x = .x, n = n) - dplyr::lag(x = .x, n = n)) / (2*n/samplingRate)
      .o <- .o + 1
    }
    return(.x)
  }
)


#' @param suffix A character string. If \code{x} is a data frame, the names of the
#'   returned data frame are equal to the names of \code{x}, each suffixed by
#'   \code{suffix}. Default is \code{""}.
#'
#' @usage \S4method{CentralDifference}{data.frame}(x, samplingRate, n, order, suffix)
#'
#' @name CentralDifference,data.frame-method
#' @rdname CentralDifference
methods::setMethod(
  f = "CentralDifference",
  signature = c(x = "data.frame"),
  definition = function(x, samplingRate = 1, n = 1, order = 1, suffix = "") {
    .diffed <- purrr::map_df(x, CentralDifference,
                             samplingRate = samplingRate, n = n, order = order)
    .named <- purrr::set_names(.diffed, stringr::str_c(names(.diffed), suffix))
    return(.named)
  }
)

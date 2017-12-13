
#' ag501: Carstens AG501 Electromagnetic Articulograph
#'
#' An R package for reading and analyzing kinematic data recorded from a Carstens
#' AG501 electromagnetic articulograph (EMA). Metadata and positional data
#' (x-, y-, and z-dimensions, and rotation) recorded from EMA sensors are
#' structured as a nested data frame. Functions are provided for filtering
#' positional data (e.g., low-pass Butterworth filter) or for calculating derived
#' time-series (e.g., velocity or acceleration profiles from either forward or
#' central differencing).
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @section Imports:
#'   \itemize{
#'     \item \code{dplyr}
#'     \item \code{lubridate}
#'     \item \code{magrittr}
#'     \item \code{methods}
#'     \item \code{purrr}
#'     \item \code{readr}
#'     \item \code{rlang}
#'     \item \code{signal}
#'     \item \code{stringi}
#'     \item \code{stringr}
#'     \item \code{tibble}
#'   }
#'
#' @section Functions:
#'   \itemize{
#'     \item \code{\link{FormatSensors}}
#'     \item \code{\link{FormatChannels}}
#'     \item \code{\link{ReadSweepMetadata}}
#'     \item \code{\link{ReadSweepData}}
#'     \item \code{\link{ReadSweep}}
#'     \item \code{\link{Butterworth}}
#'     \item \code{\link{TimeSlice}}
#'     \item \code{\link{CentralDifference}}
#'     \item \code{\link{ForwardDifference}}
#'   }
#'
#' @docType package
#' @name ag501
NULL


#' ag501: Carstens AG501 Electromagnetic Articulograph
#'
#' Analyze kinematic data recorded from a Carstens AG501
#' electromagnetic articulograph (EMA). Metadata and positional data recorded
#' from EMA sensors are structured as a nested data frame. Functions provided
#' for filtering data (e.g., low-pass Butterworth filter) or for calculating
#' derived time-series from the position data (e.g., velocity or acceleration
#' from forward or central differencing) are defined with data-first argument
#' lists, making them compatible with the tidyverse.
#'
#' @importFrom magrittr %>%
#'
#' @section Imports:
#'   \itemize{
#'     \item \code{dplyr}
#'     \item \code{lubridate}
#'     \item \code{magrittr}
#'     \item \code{methods}
#'     \item \code{purrr}
#'     \item \code{readr}
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
#'     \item \code{\link{ReadBatch}}
#'     \item \code{\link{Butterworth}}
#'     \item \code{\link{CentralDifference}}
#'     \item \code{\link{ForwardDifference}}
#'   }
#'
#' @docType package
#' @name ag501
NULL
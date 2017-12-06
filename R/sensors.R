


#' Format Sensor Names
#'
#' Format a character vector whose elements denote sensor names.
#'
#' @param sensors A character vector, names for sensors that were attached to
#'   the participant and connected to sockets in the Carstens Sensin box. The
#'   position of a name in this vector is interpreted as the number of the
#'   socket to which the sensor was connected. Default is \code{character()}.
#' @param n An atomic numeric, the number of sensors/sockets for which the
#'   Carstens Sweepsaver software recorded data. (This number is recorded in
#'   the \code{NumberOfChannels} field in the header of a \code{.pos} file.)
#'   Default is \code{16}.
#'
#' @return A character vector, derived by padding \code{sensors} with empty
#'   strings (\code{""}) to length \code{n}, and then replacing the empty
#'   strings with \code{"CH##"}, where \code{##} denotes the position of the
#'   string.
#'
#' @examples
#' FormatSensors(n = 8)
#' FormatSensors(sensors = c("HL", "HR", "TB", "TD", "TT"), n = 16)
#' FormatSensors(sensors = c("", "", "TB", "TD", "TT"), n = 16)
#'
#' @export
FormatSensors <- function(sensors = character(), n = 16) {
  .channels <- paste0("CH", sprintf("%02d", 1:n))
  .sensors <- ifelse(1:n <= length(sensors) & sensors[1:n] != "", sensors, .channels)
  return(.sensors)
}


#' Format Channel Names
#'
#' Expand a character vector whose elements denote sensor names into a character
#' vector whose elements denote names for the channels of data that were
#' recorded for those sensors.
#'
#' For each sensor, seven channels of data are recorded:
#' \itemize{
#'   \item \code{x}: position of the sensor along the antero-posterior axis
#'                   (i.e., orthogonal to the coronal plane)
#'   \item \code{y}: position of the sensor along the lateral axis
#'                   (i.e., orthogonal to the sagittal plane)
#'   \item \code{z}: position of the sensor along the vertical axis
#'                   (i.e., orthogonal to the transverse plane)
#'   \item \code{phi}: rotation angle between the \code{x}-axis and the
#'                     projection of the sensor axis onto the transverse plane
#'                     (i.e., the azimuthal angle)
#'   \item \code{theta}: rotation angle between the \code{z}-axis and the
#'                       sensor axis (i.e., the polar angle)
#'   \item \code{rms}: the root-mean-squared error between the measured
#'                     amplitudes from the articulograph's transmitter coils
#'                     and the expected amplitudes given the positional
#'                     and rotational coordinates of the sensor
#'   \item \code{extra}: an empty channel that serves as a placeholder for
#'                       future development
#' }
#' The \code{x}, \code{y}, and \code{z} axes define a three-dimensional
#' right-handed Cartesian coordinate system. The \code{phi} and \code{theta}
#' angles define the sensor axis in a spherical coordinate system.
#'
#' Because seven channels of data are recorded for each sensor, the returned
#' character vector has \code{7n} elements.
#'
#' Note that the MATLAB toolbox SMASH names the positional axes differently:
#' the antero-posterior axis is named \code{z}; the lateral axis, \code{x}; and
#' the vertical axis, \code{y}.
#'
#' @param sensors A character vector, names for sensors that were attached to
#'   the participant and connected to sockets in the Carstens Sensin box. The
#'   position of a name in this vector is interpreted as the number of the
#'   socket to which the sensor was connected. Default is \code{character()}.
#' @param n An atomic numeric, the number of sensors/sockets for which the
#'   Carstens Sweepsaver software recorded data. (This number is recorded in
#'   the \code{NumberOfChannels} field in the header of a \code{.pos} file.)
#'   Default is \code{16}.
#'
#' @return A character vector, the result of first formatting \code{sensors}
#'   (by calling \code{FormatSensors}) and then suffixing each formatted
#'   sensor name with \code{x}, \code{y}, \code{z}, \code{phi}, \code{theta},
#'   \code{rms}, and \code{extra}.
#'
#' @examples
#' FormatChannels(n = 8)
#' FormatChannels(sensors = c("HL", "HR", "TB", "TD", "TT"), n = 16)
#' FormatChannels(sensors = c("", "", "TB", "TD", "TT"), n = 16)
#'
#' @seealso \code{\link{FormatSensors}}
#'
#' @export
FormatChannels <- function(sensors = character(), n = 16) {
  FormatSensors(sensors = sensors, n = n) %>%
    purrr::map(function(.s) {paste0(.s, c("x", "y", "z", "phi", "theta", "rms", "extra"))}) %>%
    purrr::reduce(c)
}



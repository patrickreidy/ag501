

#' @return
#' \itemize{
#'   \item \code{Sweep}: A character, the name of the sweep, which identifies
#'     files with \code{.pos} and \code{.txt} extensions.
#'   \item \code{Format}: A character, the data format defined by Carstens.
#'   \item \code{SamplingRate}: A numeric, the number of samples recorded
#'     per second by the AG501.
#'   \item \code{SweepsaverBuild}: A character, the version and revision number
#'     of the Sweepsaver application, which records the amplitude values
#'     detected by the sensor coils in the AG501.
#'   \item \code{SweepsaverVersion}: A numeric, the version number of the
#'     Sweepsaver application.
#'   \item \code{SweepsaverRevision}: A numeric, the revision number of the
#'     Sweepsaver application.
#'   \item \code{SweepsaverTimestamp}: A datetime, the date and time when the
#'     sweep was recorded by the Sweepsaver application, in
#'     \code{YYYY-MM-DD HH:MM:SS} format.
#'   \item \code{CalcposBuild}: A character, the version and revision number
#'     of the Calcpos application, which transforms the amplitude values
#'     recorded by Sweepsaver into raw position data that has not been
#'     corrected for head movement.
#'   \item \code{CalcposVersion}: A numeric, the version number of the Calcpos
#'     application.
#'   \item \code{CalcposRevision}: A numeric, the revision number of the Calcpos
#'     application.
#'   \item \code{CalcposTimestamp}: A datetime, the date and time when the
#'     raw position data was computed by Calcpos.
#'   \item \code{NormposBuild}: A character, the version and revision number
#'     of the Normpos application, which transforms the raw position data by
#'     removing head movement. If the Carstens biteplane was used to compute
#'     a rotation matrix, then Normpos also rotates and translates the data
#'     so that the participant's biteplane lies within the x-y plane and so
#'     that the participant's upper incisors lie at the origin of the coordinate
#'     system.
#'   \item \code{NormposVersion}: A numeric, the version number of the Normpos
#'     application.
#'   \item \code{NormposRevision}: A numeric, the revision number of the Normpos
#'     application.
#'   \item \code{NormposTimestamp}: A datetime, the date and time when the
#'     head-correction and biteplane-rotation were applied by Normpos.
#'   \item \code{TaxonomicDistanceMean}: A numeric, the mean of the taxonomic
#'     distances, averaging across samples. The taxonomic distance of each
#'     sample is, roughly, the difference between the measured position and the
#'     expected position relative to a reference object.
#'   \item \code{TaxonomicDistanceStdDev}: A numeric, the standard deviation of
#'     taxonomic distances for the sweep.
#'   \item \code{PostFilters}: A data frame that records the low-pass Kaiser
#'     filters that were applied to the data channels of the active sensors,
#'     during postprocessing by the Calcpos and Normpos applications. Includes
#'     the following columns:
#'     \itemize{
#'       \item \code{Socket}: A numeric, the number of the socket in the Sensin
#'         box to which the physical sensor was connected.
#'       \item \code{Sensor}: A character, the name of the sensor whose data
#'         channels were low-pass filtered. See: \code{\link{FormatChannels}}.
#'       \item \code{Name}: A character, the name of the filter applied to the
#'         sensor's data channels. See: \code{\link{CarstensFilters}}.
#'       \item \code{Application}: A character, the application that applied
#'         the low-pass filter to the data.
#'     }
#' }



#' Greek Clusters
#'
#' Metadata for productions of Greek /sk/ and /sp/ clusters.
#'
#' @param includePaths If \code{TRUE}, then the paths to the \code{.pos} and
#'   \code{.txt} files in the installed package are found, and included as
#'   variables named \code{POS} and \code{TXT}, respectively.
#'
#'   If \code{FALSE}, then the paths to these files for the sweep are not
#'   included.
#'
#' @return If \code{includePaths = TRUE}, a data table with 4 rows and
#'   9 variables:
#'   \itemize{
#'     \item \code{Sweep}: the sweep number in which the cluster was produced
#'       and recorded.
#'     \item \code{Word}, \code{Cluster}, \code{Vowel}, \code{Anchor}:
#'       the WorldBet transcriptions of the target word, the word-initial
#'       fricative-stop cluster, the post-cluster vowel, and the post-vocalic
#'       singleton consonant, respectively.
#'     \item \code{Onset}, \code{Offset}:
#'       the approximate times of the onset and offset of the production
#'       of the word in the sweep. (Note: These times do not consistently
#'       correspond to any well-defined acoustic or kinematic events.)
#'     \item \code{POS}, \code{TXT}:
#'       the paths to the \code{.pos} and \code{.txt} files in the installed
#'       package for the sweeps.
#'   }
#'
#'   If \code{includePaths = FALSE}, then returned data table omits the
#'   \code{POS} and \code{TXT} variables
#'
#' @export
GreekClusters <- function(includePaths = TRUE) {
  .greek_clusters <- greek_clusters
  if (includePaths) {
    .greek_clusters <-
      .greek_clusters %>%
      dplyr::mutate(POS = purrr::map_chr(.data$Sweep, function(.x) {
        system.file(file.path("extdata", sprintf("%s.pos", .x)),
                    package = "ag501")
      })) %>%
      dplyr::mutate(TXT = purrr::map_chr(.data$Sweep, function(.x) {
        system.file(file.path("extdata", sprintf("%s.txt", .x)),
                    package = "ag501")
      }))
  }
  return(.greek_clusters)
}

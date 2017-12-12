

#' List Sweeps
#'
#' List the \code{.pos} and \code{.txt} files for sweeps. Vectorized over
#' \code{path} and \code{pattern} arguments.
#'
#' @param path A character vector of paths to directories that contain sweeps
#'   (pairs of \code{.pos} and \code{.txt} files that share the same basename).
#'   Tilde expansion (see \code{\link{path.expand}}) is performed.
#'   Default value corresponds to the present working directory.
#' @param pattern An optional vector of \code{\link[=regex]{regular expressions}}.
#'   The length of this vector should be equal to the length of the \code{path}
#'   argument, or equal to \code{1}, in which case the same regular expression
#'   is used to match sweeps in each directory listed in \code{path}.
#'   Only sweeps whose basenames match the regular expression(s) will be
#'   returned. Default value matches all available sweeps. Patterns are always
#'   case-sensitive.
#' @param fullNames If \code{TRUE}, then the full paths for the \code{.pos}
#'   and \code{.txt} files are returned. If \code{FALSE}, then the basenames
#'   for these files are returned. Default is \code{TRUE}.
#' @param simplify If \code{TRUE}, then the table of sweeps found for each element
#'   in \code{path} are row-binded together into a single data table.
#'   If \code{FALSE}, then the tables of sweeps found for the elements in \code{path}
#'   are returned as elements of a list.
#' @param recursive If \code{TRUE}, then the search recurses into subdirectories of
#'   elements of \code{param}. If \code{FALSE}, then the search only looks
#'   immediately in the directories listed in \code{path}.
#' @param dropMissing If \code{TRUE}, then sweeps that are missing either a
#'   \code{.pos} or a \code{.txt} file are not listed. If \code{FALSE}, then
#'   sweeps with a missing file are listed.
#'
#' @return If \code{simplify = TRUE}, a data table with the following variables:
#'   \itemize{
#'     \item \code{Sweep}: the basenames of the sweeps.
#'     \item \code{POS}: the \code{.pos} files for the sweeps.
#'     \item \code{TXT}: the \code{.txt} files for the sweeps.
#'   }
#'
#'   If \code{simplify = FALSE}, a list with as many elements as \code{path}.
#'   Each element is a data table as described immediately above.
#'
#' @export
ListSweeps <- function(path = ".", pattern = "*", fullNames = TRUE, simplify = TRUE, recursive = FALSE, dropMissing = TRUE) {
  .ListSweepsInSinglePath <- function(.path, .pattern, .fn, .rec) {
    .files <- list.files(path = .path, pattern = .pattern, full.names = .fn, recursive = .rec)
    .pos <- grep(pattern = "\\.pos$", x = .files, value = TRUE)
    .pos_tbl <- tibble::tibble(Sweep = stringr::str_replace(basename(.pos), "\\.pos$", ""),
                               POS = .pos)
    .txt <- grep(pattern = "\\.txt$", x = .files, value = TRUE)
    .txt_tbl <- tibble::tibble(Sweep = stringr::str_replace(basename(.txt), "\\.txt", ""),
                               TXT = .txt)
    .sweeps <-
      dplyr::full_join(x = .pos_tbl, y = .txt_tbl, by = "Sweep") %>%
      dplyr::mutate(N = ifelse(!is.na(.data$POS), 1, 0) + ifelse(!is.na(.data$TXT), 1, 0)) %>%
      dplyr::arrange(dplyr::desc(.data$N), .data$Sweep) %>%
      dplyr::select(-.data$N)
    if (dropMissing) {
      .sweeps <- dplyr::filter(.sweeps, !is.na(POS), !is.na(TXT))
    }
    .listed <- purrr::set_names(list(.sweeps), .path)
    return(.listed)
  }
  if (length(pattern) != length(path) & length(pattern) != 1) {
    .stop_msg <- sprintf("length of @pattern (%d) must be either length of @path (%d) or 1.",
                         length(pattern), length(path))
    stop(.stop_msg)
  } else if (length(pattern) == 1) {
    pattern <- rep(pattern, length(path))
  }
  .sweeps <- purrr::map2(path, pattern, .ListSweepsInSinglePath,
                         .fn = fullNames, .rec = recursive)
  if (simplify) {
    if (length(.sweeps) == 1) {
      .sweeps <- purrr::flatten_df(.sweeps)
    } else {
      .sweeps <- purrr::reduce(.sweeps, dplyr::bind_rows)
    }
  } else {
    .sweeps <- purrr::flatten(.sweeps)
  }
  return(.sweeps)
}

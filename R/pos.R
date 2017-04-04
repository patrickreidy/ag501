

parseDbl <- function(header, field, name) {
  header %>%
    stringr::str_subset(pattern = field) %>%
    stringr::str_split(pattern = "=") %>%
    purrr::map_chr(2) %>%
    as.numeric() %>%
    tibble::as_tibble() %>%
    purrr::set_names(name)
}


parseBuild <- function(header, field, prefix) {
  header %>%
    stringr::str_subset(pattern = field) %>%
    stringr::str_split(pattern = "=") %>%
    purrr::map_chr(2) %>%
    tibble::as_tibble() %>%
    dplyr::rename(Build = value) %>%
    dplyr::mutate(
      Version = Build %>%
        stringr::str_extract(pattern = "[0-9]+\\.[0-9]+") %>%
        as.numeric()
    ) %>%
    dplyr::mutate(
      Revision = Build %>%
        stringr::str_extract(pattern = "[0-9]+$") %>%
        as.numeric()
    ) %>%
    purrr::set_names(paste(prefix, names(.), sep = ""))
}


parseTimestamp <- function(header, field, name) {
  header %>%
    stringr::str_subset(pattern = field) %>%
    stringr::str_split(pattern = "=") %>%
    purrr::map_chr(2) %>%
    lubridate::ymd_hms(tz = Sys.timezone()) %>%
    tibble::as_tibble() %>%
    purrr::set_names(name)
}


parseTaxonomic <- function(header, statistic) {
  .normpos_version <- parseBuild(header, "normpos.version", "") %>%
    dplyr::select(Version) %>%
    as.numeric()
  .field <- stringr::str_interp("Taxonomic_Distance_${statistic}")
  .taxonomic <- ifelse(
    .normpos_version >= 2.4,
    parseDbl(header, .field, "x") %>% as.numeric(),
    as.numeric(NA)
  )
  tibble::tibble(.taxonomic) %>%
    purrr::set_names(
      .field %>% stringr::str_replace_all(pattern = "_", replacement = "")
    )
}


parseCalcposFilter <- function(filter) {
  tibble::tibble(
    Name = filter %>%
      stringr::str_split(pattern = "=") %>%
      purrr::map_chr(2),
    Application = "Calcpos"
  )
}


parseNormposFilter <- function(filter, sensors) {
  tibble::tibble(
    Socket = filter %>%
      stringr::str_split(pattern = "=") %>%
      purrr::map_chr(2) %>%
      stringr::str_split(pattern = ",", simplify = TRUE) %>%
      as.numeric(),
    Sensor = sensors[Socket],
    Name = filter %>%
      stringr::str_split(pattern = "\\.") %>%
      purrr::map_chr(2) %>%
      stringr::str_split(pattern = "=") %>%
      purrr::map_chr(1),
    Application = "Normpos"
  )
}





#' Read the Header of the .pos File for a Sweep
#'
#' Read the header of a \code{.pos} file, which contains the metadata for the
#' positional data recorded during a sweep.
#'
#' @param file A character string, the path to a \code{.pos} file.
#' @param sensors A character vector, names for sensors that were attached to
#'   the participant and connected to sockets in the Carstens Sensin box. The
#'   position of a name in this vector is interpreted as the number of the
#'   socket to which the sensor was connected. Default is \code{c("")}.
#'
#' @return A data tibble with the following variables:
#' \itemize{
#'   \item \code{Sweep}
#'   \item \code{Format}
#'   \item \code{SamplingRate}
#'   \item \code{SweepsaverBuild}
#'   \item \code{SweepsaverVersion}
#'   \item \code{SweepsaverRevision}
#'   \item \code{SweepsaverTimestamp}
#'   \item \code{CalcposBuild}
#'   \item \code{CalcposVersion}
#'   \item \code{CalcposRevision}
#'   \item \code{CalcposTimestamp}
#'   \item \code{NormposBuild}
#'   \item \code{NormposVersion}
#'   \item \code{NormposRevision}
#'   \item \code{NormposTimestamp}
#'   \item \code{TaxonomicDistanceMean}
#'   \item \code{TaxonomicDistanceStdDev}
#'   \item \code{PostFilters}
#' }
#'
#' @export
ReadSweepMetadata <- function(file, sensors = c("")) {
  # Read the lines of the `file`. Select just the header.
  .header <- readr::read_lines(file) %>%
    .[1:(match("", .) - 1)]
  # Format the names of the `sensors`.
  .sensors <- sensors %>%
    FormatSensors(n = parseDbl(.header, "NumberOfChannels", "x") %>%
                    as.numeric())
  # Parse the filters applied by cs5normpos.
  .normpos_filters <- .header %>%
    stringr::str_subset(pattern = "normpos") %>%
    stringi::stri_subset_regex(pattern = "version|timestamp|Taxonomic_Distance",
                               negate = TRUE) %>%
    purrr::map(parseNormposFilter, sensors = .sensors) %>%
    purrr::reduce(dplyr::bind_rows)
  # Parse the filter applied by cs5calcpos.
  .calcpos_filters <- .header %>%
    stringr::str_subset(pattern = "calcpos") %>%
    stringi::stri_subset_regex(pattern = "version|timestamp", negate = TRUE) %>%
    parseCalcposFilter() %>%
    dplyr::sample_n(size = nrow(.normpos_filters), replace = TRUE) %>%
    dplyr::bind_cols(dplyr::select(.normpos_filters, Socket, Sensor)) %>%
    dplyr::select(Socket, Sensor, Name, Application)
  # Bind the normpos-filters and the calcpos-filters.
  .post_filters <- .calcpos_filters %>%
    dplyr::bind_rows(.normpos_filters) %>%
    dplyr::arrange(Socket, Application)
  .metadata <- tibble::tibble(
    Sweep = basename(file) %>%
      stringr::str_replace(pattern = ".pos", replacement = ""),
    Format = .header %>%
      stringr::str_subset(pattern = "AG50xDATA")
  ) %>%
    dplyr::bind_cols(
      parseDbl(header = .header,
               field = "SamplingFrequency",
               name = "SamplingRate"),
      parseBuild(header = .header,
                 field = "sweepsaver.version",
                 prefix = "Sweepsaver"),
      parseTimestamp(header = .header,
                     field = "recorded",
                     name = "SweepsaverTimestamp"),
      parseBuild(header = .header,
                 field = "calcpos.version",
                 prefix = "Calcpos"),
      parseTimestamp(header = .header,
                     field = "calcpos.timestamp",
                     name = "CalcposTimestamp"),
      parseBuild(header = .header,
                 field = "normpos.version",
                 prefix = "Normpos"),
      parseTimestamp(header = .header,
                     field = "normpos.timestamp",
                     name = "NormposTimestamp"),
      parseTaxonomic(header = .header,
                     statistic = "Mean"),
      parseTaxonomic(header = .header,
                     statistic = "StdDev")
    ) %>%
    dplyr::mutate(PostFilters = list(.post_filters))
  return(.metadata)
}





#' Read the Contents of the .txt File for a Sweep
#'
#' Read the contents of a \code{.txt} file, which contains the positional data
#' recorded during a sweep.
#'
#' @param file A character string, the path to a \code{.pos} file.
#' @param sensors A character vector, names for sensors that were attached to
#'   the participant and connected to sockets in the Carstens Sensin box. The
#'   position of a name in this vector is interpreted as the number of the
#'   socket to which the sensor was connected. Default is \code{c("")}.
#'
#' @return A data tibble, the names of whose variables are determined by
#'   expanding the \code{sensors} argument through \code{FormatChannels} and
#'   then dropping the channels suffixed by \code{extra} (i.e., the placeholder
#'   channels).
#'
#' @seealso \code{\link{FormatChannels}}
#'
#' @export
ReadSweepData <- function(file, sensors = c("")) {
  suppressMessages(readr::read_delim(file, delim = "\t", col_names = FALSE)) %>%
    purrr::set_names(FormatChannels(sensors = sensors, n = ncol(.)/7)) %>%
    dplyr::select(-dplyr::ends_with("extra"))
}





#' Read the Positional Data Recorded During a Sweep
#'
#' Read the head-corrected positional data from a \code{.txt} file and its
#' metadata from the header of a \code{.pos} file.
#'
#' @param sweep A character string, the basename of associated \code{.txt} and
#'   \code{.pos} files that are generated from post-processing a sweep. That is,
#'   there should be files \code{<sweep>.txt} and \code{<sweep>.pos} located in
#'   \code{path}.
#' @param path A character string, the path to the directory containing the
#'   \code{.txt} and \code{.pos} files.
#' @param sensors A character vector, names for sensors that were attached to
#'   the participant and connected to sockets in the Carstens Sensin box. The
#'   position of a name in this vector is interpreted as the number of the
#'   socket to which the sensor was connected. Default is \code{c("")}.
#'
#' @return A data tibble, with the following variables:
#' \itemize{
#'   \item all variables returned by \code{ReadSweepMetadata}
#'   \item \code{TimeData}: a list of data tibbles, the times at which position
#'         data samples were recorded
#'   \item \code{PositionData}: a list of data tibbles, the position data
#'         returned by \code{ReadSweepData} for the active sensors that were
#'         recorded during the sweep and then post-processed with Calcpos and
#'         Normpos
#' }
#'
#' @seealso \code{\link{ReadSweepMetadata}}, \code{\link{ReadSweepData}}
#'
#' @export
ReadSweep <- function(sweep, path, sensors = c("")) {
  .metadata <- list.files(path = path,
                          pattern = stringr::str_interp("${sweep}.pos"),
                          full.names = TRUE) %>%
    ReadSweepMetadata(sensors = sensors)
  .data <- list.files(path = path,
                      pattern = stringr::str_interp("${sweep}.txt"),
                      full.names = TRUE) %>%
    ReadSweepData(sensors = sensors)
  .active_cols <- .metadata %>%
    dplyr::select(PostFilters) %>%
    purrr::map(1) %>%
    purrr::map("Sensor") %>%
    purrr::reduce(c) %>%
    unique() %>%
    purrr::map(dplyr::starts_with, vars = names(.data)) %>%
    purrr::reduce(c)
  .pos_data <- .metadata %>%
    dplyr::mutate(PositionData = list(.data %>% dplyr::select(.active_cols))) %>%
    dplyr::mutate(TimeData = purrr::map2(
      SamplingRate, PositionData,
      ~ tibble::tibble(Time = seq(from = 0, by = 1/.x, length.out = nrow(.y)))
    )) %>%
    dplyr::select(-dplyr::one_of(c("TimeData", "PositionData")), TimeData, PositionData)
  return(.pos_data)
}




#' Read the Positional Data for a Batch of Sweeps
#'
#' Read head-corrected positional data from \code{.txt} files and metadata
#' from the headers of \code{.pos} files.
#'
#' @param path A character string, the path to the directory containing the
#'   \code{.txt} and \code{.pos} files.
#' @param sensors A character vector, names for sensors that were attached to
#'   the participant and connected to sockets in the Carstens Sensin box. The
#'   position of a name in this vector is interpreted as the number of the
#'   socket to which the sensor was connected. Default is \code{c("")}.
#'
#' @return A data tibble with the same variables as those returned by
#'   \code{ReadPos}.
#'
#' @seealso \code{\link{ReadSweep}}
#'
#' @export
ReadBatch <- function(path, sensors = c("")) {
  .pos <- list.files(path = path, pattern = "\\.pos$") %>%
    stringr::str_replace(pattern = "\\.pos$", replacement = "")
  .txt <- list.files(path = path, pattern = "\\.txt$") %>%
    stringr::str_replace(pattern = "\\.txt$", replacement = "")
  intersect(.pos, .txt) %>%
    purrr::map(ReadSweep, path = path, sensors = sensors) %>%
    purrr::reduce(dplyr::bind_rows)
}



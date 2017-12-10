

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
  .parsed <- header %>%
    stringr::str_subset(pattern = field) %>%
    stringr::str_split(pattern = "=") %>%
    purrr::map_chr(2) %>%
    tibble::as_tibble() %>%
    dplyr::rename(Build = .data$value) %>%
    dplyr::mutate(
      Version = stringr::str_extract(string = .data$Build,
                                     pattern = "[0-9]+\\.[0-9]+"),
      Version = as.numeric(.data$Version)
    ) %>%
    dplyr::mutate(
      Revision = stringr::str_extract(string = .data$Build,
                                      pattern = "[0-9]+$"),
      Revision = as.numeric(.data$Revision)
    )
  .parsed_named <- purrr::set_names(.parsed, nm = stringr::str_c(prefix, names(.parsed)))
  return(.parsed_named)
}


parseTimestamp <- function(header, field, name) {
  .tz <- Sys.timezone()
  if (! (.tz %in% OlsonNames())) {
    Sys.setenv(TZ = "UTC")
  }
  .timezone <-
    header %>%
    stringr::str_subset(pattern = field) %>%
    stringr::str_split(pattern = "=") %>%
    purrr::map_chr(2) %>%
    lubridate::ymd_hms(tz = Sys.timezone()) %>%
    tibble::as_tibble() %>%
    purrr::set_names(name)
  if (! (.tz %in% OlsonNames())) {
    Sys.unsetenv("TZ")
  }
  return(.timezone)
}


parseTaxonomic <- function(header, statistic) {
  .normpos_version <- parseBuild(header, "normpos.version", "") %>%
    dplyr::select(.data$Version) %>%
    as.numeric()
  .field <- stringr::str_interp("Taxonomic_Distance_${statistic}")
  .taxonomic <- ifelse(.normpos_version >= 2.4,
                       as.numeric(parseDbl(header, .field, "x")),
                       as.numeric(NA))
  .tax_tibble <- purrr::set_names(tibble::tibble(.taxonomic),
                                    stringr::str_replace_all(.field, "_", ""))
  return(.tax_tibble)
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
    Sensor = sensors[.data$Socket],
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
#' @param pos A character string, the path to a \code{.pos} file.
#' @param sensors A character vector, names for sensors that were attached to
#'   the participant and connected to sockets in the Carstens Sensin box. The
#'   position of a name in this vector is interpreted as the number of the
#'   socket to which the sensor was connected. Default is \code{character()}.
#'   See \code{\link{FormatSensors}}.
#'
#' @return A data tibble with the following variables:
#' @template sweep-metadata
#'
#' @export
ReadSweepMetadata <- function(pos, sensors = character()) {
  # Read the lines of the `pos`. Select just the header.
  .lines <- readr::read_lines(pos)
  .header <- .lines[1:(match("", .lines) - 1)]
  # Format the names of the `sensors`.
  .n <- as.numeric(parseDbl(.header, field = "NumberOfChannels", name = "x"))
  .sensors <- FormatSensors(sensors = sensors, n = .n)
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
    dplyr::bind_cols(dplyr::select(.normpos_filters, .data$Socket, .data$Sensor)) %>%
    dplyr::select(.data$Socket, .data$Sensor, .data$Name, .data$Application)
  # Bind the normpos-filters and the calcpos-filters.
  .post_filters <- .calcpos_filters %>%
    dplyr::bind_rows(.normpos_filters) %>%
    dplyr::arrange(.data$Socket, .data$Application)
  .metadata <- tibble::tibble(
    Sweep = basename(pos) %>%
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
#' @param txt A character string, the path to a \code{.txt} file.
#' @param sensors A character vector, names for sensors that were attached to
#'   the participant and connected to sockets in the Carstens Sensin box. The
#'   position of a name in this vector is interpreted as the number of the
#'   socket to which the sensor was connected. Default is \code{character()}.
#'   See \code{\link{FormatSensors}}.
#' @param n A numeric, the number of data channels that are recorded
#'   per sensor. Default is \code{7}, which is consistent with Carstens data
#'   format \code{AG50xDATA_V003}, under which each sensor comprises data
#'   channels: \code{x}, \code{y}, \code{z}, \code{phi}, \code{theta},
#'   \code{rms}, \code{extra}.
#'   See \code{\link{FormatChannels}}.
#'
#' @return A data tibble, the names of whose variables are determined by
#'   expanding the \code{sensors} argument through \code{\link{FormatChannels}}
#'   and dropping the channels suffixed by \code{extra}, which contain no
#'   meaningful data.
#'
#' @export
ReadSweepData <- function(txt, sensors = character(), n = 7) {
  .channel_data <- suppressMessages(readr::read_delim(txt, delim = "\t", col_names = FALSE))
  .channel_names <- FormatChannels(sensors = sensors, n = ncol(.channel_data)/n)
  .sweep_data <- purrr::set_names(.channel_data, .channel_names) %>%
    dplyr::select(-dplyr::ends_with("extra"))
  return(.sweep_data)
}





#' Read the Positional Data Recorded During a Sweep
#'
#' Read the head-corrected (and biteplane-rotated) positional data from a
#' \code{.txt} file and its metadata from the header of a \code{.pos} file.
#'
#' @param pos A character string, the path to a \code{.pos} file.
#'   Default is \code{""}.
#' @param txt A character string, the path to a \code{.txt} file.
#'   Default is \code{""}.
#' @param sweep A character string, the path to and basename (without extension)
#'   of a sweep. \code{stringr::str_interp("${sweep}.pos")} should point to
#'   a \code{.pos} file whose header contains the metadata about the sweep.
#'   \code{stringr::str_interp("${sweep}.txt")} should point to a \code{.txt}
#'   file that contains the positional data recorded for the sweep.
#'   Default is \code{NULL}. If not \code{NULL}, then \code{-.pos} and
#'   \code{-.txt} completions of \code{sweep} override any provided values of
#'   the \code{pos} and \code{txt} arguments.
#' @param sensors A character vector, names for sensors that were attached to
#'   the participant and connected to sockets in the Carstens Sensin box. The
#'   position of a name in this vector is interpreted as the number of the
#'   socket to which the sensor was connected. Default is \code{character()}.
#'   See \code{\link{FormatSensors}} and \code{\link{FormatChannels}}.
#'
#' @return A nested data frame with the following variables:
#' @template sweep-metadata
#' @template sweep-data
#'
#' @export
ReadSweep <- function(pos = "", txt = "", sweep = NULL, sensors = character()) {
  if (!is.null(sweep)) {
    pos <- stringr::str_interp("${sweep}.pos")
    txt <- stringr::str_interp("${sweep}.txt")
  }
  .metadata <- ReadSweepMetadata(pos = pos, sensors = sensors)
  .active_chs <-
    .metadata %>%
    dplyr::pull(.data$PostFilters) %>%
    purrr::flatten_df() %>%
    dplyr::filter(.data$Application == "Normpos") %>%
    dplyr::pull(.data$Sensor) %>%
    paste(collapse = "|")
  .channel_data <-
    ReadSweepData(txt = txt, sensors = sensors) %>%
    dplyr::select(dplyr::matches(.active_chs))
  .time_samples <- seq(from = 0, by = 1/dplyr::pull(.metadata, .data$SamplingRate),
                       length.out = nrow(.channel_data))
  .time_data <- tibble::tibble(Time = .time_samples)
  .position_data <- dplyr::bind_cols(.time_data, .channel_data)
  .sweep <- dplyr::mutate(.metadata,
                          PositionData = list(.position_data))
  return(.sweep)
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
#' @return A nested data frame, the result of calling \code{\link{ReadSweep}}
#'   on each sweep that has a \code{.pos} and a \code{.txt} file in \code{path},
#'   and then row-binding them together. The returned data frame has the
#'   following columns:
#' @template sweep-metadata
#' @template sweep-data
#'
#' @noRd
ReadBatch <- function(path, sensors = c("")) {
  .pos_files <-
    path %>%
    list.files(pattern = "\\.pos$") %>%
    stringr::str_replace(pattern = "\\.pos$", replacement = "")
  .txt_files <-
    path %>%
    list.files(pattern = "\\.txt$") %>%
    stringr::str_replace(pattern = "\\.txt$", replacement = "")
  .sweeps <-
    path %>%
    file.path(dplyr::intersect(.pos_files, .txt_files)) %>%
    purrr::map(ReadSweep, sensors = sensors) %>%
    purrr::reduce(dplyr::bind_rows)
  return(.sweeps)
}



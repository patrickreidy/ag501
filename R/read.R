

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
#' \code{ReadSweepMetadata} is vectorized over both \code{pos} and \code{sensors}.
#' If \code{pos} is a vector of more than one path to a \code{.pos} file, and
#' if \code{sensors} is a character vector, then the elements of \code{sensors}
#' are recycled as the sensor-names for each \code{.pos} file. To specify
#' different sensor names for each \code{.pos} file, \code{sensors} should be
#' a list with the same length as \code{pos} and whose elements are character
#' vectors.
#'
#' @param pos A character vector of paths to \code{.pos} files.
#' @param sensors A character vector or list of character vectors, the
#'   names for sensors that were attached to the participant and connected to
#'   sockets in the Carstens Sensin box. The position of a name in this vector
#'   is interpreted as the number of the socket to which the sensor was connected.
#'   Default is \code{character()}. See \code{\link{FormatSensors}}.
#'
#' @return If \code{pos} is a single character string, then a single-row data
#'   table with the following variables:
#'   @template sweep-metadata
#'
#'   If \code{pos} is a character vector with more than one element, then a list
#'   whose elements are data tables structured as described above.
#'
#' @export
ReadSweepMetadata <- function(pos, sensors = character()) {
  .ReadSinglePOS <- function(.pos, .sens = sensors) {
    if (!is.character(.pos)) {
      stop(".ReadSinglePOS requires a character string as the .pos argument")
    }
    if (length(.pos) > 1) {
      stop(".ReadSinglePOS can only handle a single path in the .pos argument")
    }
    if (!is.character(.sens)) {
      stop(".ReadSinglePOS requires a character vector as the .sens argument ")
    }
    # Read the lines of the `.pos` file. Select just the header.
    .lines <- readr::read_lines(.pos)
    .header <- .lines[1:(match("", .lines) - 1)]
    # Format the names of the `sensors`.
    .n <- as.numeric(parseDbl(.header, field = "NumberOfChannels", name = "x"))
    .sensors <- FormatSensors(sensors = .sens, n = .n)
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
      Sweep = basename(.pos) %>%
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
  if (length(pos) > 1) {
    if (is.character(sensors)) {
      sensors <- rep(list(sensors), times = length(pos))
    } else {
      if (is.list(sensors)) {
        if (length(sensors) == 1) {
          sensors <- rep(sensors, times = length(pos))
        } else if (length(sensors) != length(pos)) {
          .stop_msg <- sprintf("length of @sensors list (%d) must be either length of @pos (%d) or 1.",
                               length(sensors), length(pos))
          stop(.stop_msg)
        }
      }
    }
    .sweep_metadata <- purrr::map2(pos, sensors, .ReadSinglePOS)
  } else {
    .sweep_metadata <- .ReadSinglePOS(pos, sensors)
  }
  return(.sweep_metadata)
}





#' Read the Contents of the .txt File for a Sweep
#'
#' Read the contents of a \code{.txt} file, which contains the positional data
#' recorded during a sweep.
#'
#' \code{ReadSweepData} is vectorized over both \code{txt} and \code{sensors}.
#' If \code{txt} is a vector of more than one path to a \code{.txt} file, and
#' if \code{sensors} is a character vector, then the elements of \code{sensors}
#' are recycled as the sensor-names for each \code{.txt} file. To specify
#' different sensor names for each \code{.txt} file, \code{sensors} should be
#' a list with the same length as \code{txt} and whose elements are character
#' vectors.
#'
#' @param txt A character vector of paths to \code{.txt} files.
#' @param sensors A character vector or list of character vectors, the
#'   names for sensors that were attached to the participant and connected to
#'   sockets in the Carstens Sensin box. The position of a name in this vector
#'   is interpreted as the number of the socket to which the sensor was connected.
#'   Default is \code{character()}. See \code{\link{FormatSensors}}.
#' @param n A numeric, the number of data channels that are recorded
#'   per sensor. Default is \code{7}, which is consistent with Carstens data
#'   format \code{AG50xDATA_V003}, under which each sensor comprises data
#'   channels: \code{x}, \code{y}, \code{z}, \code{phi}, \code{theta},
#'   \code{rms}, \code{extra}.
#'   See \code{\link{FormatChannels}}.
#' @param dropExtra If \code{TRUE}, then the \code{extra} channel for each sensor
#'   is dropped from the data, since this channel is only a placeholder for
#'   future development by Carstens and contains no data.
#'   If \code{FALSE}, then the \code{extra} channels are kept.
#'
#' @return If \code{txt} is a single character string, then a data table,
#'   the names of whose variables are determined by expanding the \code{sensors}
#'   argument through \code{\link{FormatChannels}}, and optionally dropping
#'   the channels suffixed by \code{extra}.
#'
#'   If \code{txt} is a character vector with more than one element, then a list
#'   whose elements are data tables structured as described above.
#'
#' @export
ReadSweepData <- function(txt, sensors = character(), n = 7, dropExtra = TRUE) {
  .ReadSingleTXT <- function(.txt, .sens = sensors, .n = n, .drop = dropExtra) {
    if (!is.character(.txt)) {
      stop(".ReadSingleTXT requires a character string as the .txt argument")
    }
    if (length(.txt) > 1) {
      stop(".ReadSingleTXT can only handle a single path in the .txt argument")
    }
    if (!is.character(.sens)) {
      stop(".ReadSingleTXT requires a character vector as the .sens argument ")
    }
    .channel_data <- suppressMessages(readr::read_delim(.txt, delim = "\t", col_names = FALSE))
    .channel_names <- FormatChannels(sensors = .sens, n = ncol(.channel_data)/.n)
    .sweep_data <- purrr::set_names(.channel_data, .channel_names)
    if (dropExtra) {
      .sweep_data <- dplyr::select(.sweep_data, -dplyr::ends_with("extra"))
    }
    return(.sweep_data)
  }
  if (length(txt) > 1) {
    if (is.character(sensors)) {
      sensors <- rep(list(sensors), times = length(txt))
    } else {
      if (is.list(sensors)) {
        if (length(sensors) == 1) {
          sensors <- rep(sensors, times = length(txt))
        } else if (length(sensors) != length(txt)) {
          .stop_msg <- sprintf("length of @sensors list (%d) must be either length of @txt (%d) or 1.",
                               length(sensors), length(txt))
          stop(.stop_msg)
        }
      }
    }
    .sweep_metadata <- purrr::map2(txt, sensors, .ReadSingleTXT, .n = n, .drop = dropExtra)
  } else {
    .sweep_metadata <- .ReadSingleTXT(txt, sensors, .n = n, .drop = dropExtra)
  }
  return(.sweep_metadata)
}





#' Read the Positional Data Recorded During a Sweep
#'
#' Read the head-corrected (and biteplane-rotated) positional data from a
#' \code{.txt} file and its metadata from the header of a \code{.pos} file.
#'
#' \code{ReadSweep} is vectorized over \code{pos}, \code{txt}, and \code{sensors}.
#' \code{pos} and \code{txt} should be vectors of the same length ordered such
#' that the positions on these vectors correspond to the \code{.pos} and \code{.txt}
#' files, respectively, of the same sweep. Alternatively, multiple sweeps can
#' be specified by passing a character vector to the \code{sweep} argument.
#' If multiple sweeps are specified, and if \code{sensors} is a character vector,
#' then the elements of \code{sensors} are recycled as the sensor-names for each
#' sweep. To specify different sensor names for each sweep, \code{sensors} should be
#' a list with the same length as the number of sweeps and whose elements are
#' character vectors.
#'
#' @param pos A character vector of paths to \code{.pos} files.
#' @param pos A character vector of paths to \code{.txt} files.
#' @param sweep A character vector, the path to and basenames (without extensions)
#'   of sweeps. The names of the \code{.pos} and \code{.txt} files for the
#'   sweeps will be constructed by appending the appropriate file extensions.
#'   Default is \code{NULL}. If not \code{NULL}, then \code{-.pos} and
#'   \code{-.txt} completions of \code{sweep} override any provided values of
#'   the \code{pos} and \code{txt} arguments.
#' @param sensors A character vector, names for sensors that were attached to
#'   the participant and connected to sockets in the Carstens Sensin box. The
#'   position of a name in this vector is interpreted as the number of the
#'   socket to which the sensor was connected. Default is \code{character()}.
#'   See \code{\link{FormatSensors}} and \code{\link{FormatChannels}}.
#' @param dropExtra If \code{TRUE}, then the \code{extra} channel for each sensor
#'   is dropped from the data, since this channel is only a placeholder for
#'   future development by Carstens and contains no data.
#'   If \code{FALSE}, then the \code{extra} channels are kept.
#' @param simplify If \code{TRUE}, then the data tables for the sweeps
#'   are row-binded together into a single data table.
#'   If \code{FALSE}, then the data tables for the sweeps are returned as
#'   elements of a list.
#'
#' @return A listnested data frame with the following variables:
#' @template sweep-metadata
#' @template sweep-data
#'
#' @export
ReadSweep <- function(pos = "", txt = "", sweep = NULL, sensors = character(), dropExtra = TRUE, simplify = TRUE) {
  .ReadSingleSweep <- function(.pos, .txt, .sens = sensors, .drop = dropExtra) {
    if (!is.character(.pos)) {
      stop(".ReadSingleSweep requires a character string as the .pos argument")
    }
    if (length(.pos) > 1) {
      stop(".ReadSingleSweep can only handle a single path in the .pos argument")
    }
    if (!is.character(.txt)) {
      stop(".ReadSingleSweep requires a character string as the .txt argument")
    }
    if (length(.txt) > 1) {
      stop(".ReadSingleSweep can only handle a single path in the .txt argument")
    }
    if (!is.character(.sens)) {
      stop(".ReadSingleSweep requires a character vector as the .sens argument ")
    }
    .metadata <- ReadSweepMetadata(pos = .pos, sensors = .sens)
    .active_chs <-
      .metadata %>%
      dplyr::pull(.data$PostFilters) %>%
      purrr::flatten_df() %>%
      dplyr::filter(.data$Application == "Normpos") %>%
      dplyr::pull(.data$Sensor) %>%
      paste(collapse = "|")
    .channel_data <-
      ReadSweepData(txt = .txt, sensors = .sens, dropExtra = .drop) %>%
      dplyr::select(dplyr::matches(.active_chs))
    .time_samples <- seq(from = 0, by = 1/dplyr::pull(.metadata, .data$SamplingRate),
                         length.out = nrow(.channel_data))
    .time_data <- tibble::tibble(Time = .time_samples)
    .position_data <- dplyr::bind_cols(.time_data, .channel_data)
    .sweep <- list(dplyr::mutate(.metadata, PositionData = list(.position_data)))
    .named <- purrr::set_names(.sweep, stringr::str_replace(basename(.pos), "\\.pos$", ""))
    return(.named)
  }
  if (!is.null(sweep)) {
    pos <- sprintf(fmt = "%s.pos", sweep)
    txt <- sprintf(fmt = "%s.txt", sweep)
  }
  if (length(pos) != length(txt)) {
    .stop_msg <- sprintf("The length of @pos (%d) must be equal to the length of @txt (%d)",
                         length(pos), length(txt))
    stop(.stop_msg)
  }

  if (is.character(sensors)) {
    sensors <- list(sensors)
  }
  if (length(sensors) == 1) {
    sensors <- rep(sensors, times = length(txt))
  } else if (length(sensors) != length(txt)) {
    .stop_msg <- sprintf("length of @sensors list (%d) must be either length of @txt (%d) or 1.",
                         length(sensors), length(txt))
    stop(.stop_msg)
  }
  .sweeps <- purrr::pmap(list(pos, txt, sensors), function(.p, .t, .s) {
    .ReadSingleSweep(.pos = .p, .txt = .t, .sens = .s, .drop = dropExtra)
  })

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




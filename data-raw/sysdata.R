


# Internalize the filters that come packaged with the Carstens control server.
carstens_filters <- list.files(path = "data-raw/filters", full.names = TRUE) %>%
  tibble::as_tibble() %>%
  dplyr::rename(Path = value) %>%
  dplyr::mutate(Name = basename(Path) %>%
                  stringr::str_replace(pattern = "\\.txt$", replacement = "")) %>%
  dplyr::mutate(PassbandLimit = Name %>%
                  stringr::str_split(pattern = "_") %>%
                  purrr::map_chr(4) %>%
                  as.numeric()) %>%
  dplyr::mutate(StopbandLimit = Name %>%
                  stringr::str_split(pattern = "_") %>%
                  purrr::map_chr(5) %>%
                  as.numeric()) %>%
  dplyr::mutate(Attenuation = Name %>%
                  stringr::str_split(pattern = "_") %>%
                  purrr::map_chr(6) %>%
                  as.numeric()) %>%
  dplyr::mutate(SamplingRate = Name %>%
                  stringr::str_split(pattern = "_") %>%
                  purrr::map_chr(7) %>%
                  as.numeric()) %>%
  dplyr::mutate(DifferenceOrder = ifelse(stringr::str_detect(Name, "acc"), 2,
                                         ifelse(stringr::str_detect(Name, "vel"), 1,
                                           0))) %>%
  dplyr::mutate(Window = purrr::map(Path, ~ as.numeric(readr::read_lines(.x)))) %>%
  dplyr::select(-Path) %>%
  dplyr::arrange(DifferenceOrder, PassbandLimit, StopbandLimit, SamplingRate)










FindTagTime <- function(sweep, tag) {
  suppressMessages(
    readr::read_delim(file = sprintf(file.path("data-raw", "greek-clusters", "%s-tags.txt"), sweep),
                      delim = "\t", col_names = c("Time", "Stream", "Tag"))
  ) %>%
    dplyr::filter(Stream == "wav", Tag == sprintf("Tag%d", tag)) %>%
    dplyr::select(Time) %>%
    purrr::simplify()
}

greek_clusters <-
  tidyr::crossing(
    tibble::tibble(Sweep = c("0021", "0022")),
    tibble::tibble(Word = c("skavi", "spala"),
                   Cluster = c("sk", "sp"),
                   Vowel = "a",
                   Anchor = c("v", "l"))
  ) %>%
  dplyr::arrange(Word, Sweep) %>%
  dplyr::mutate(Production = purrr::map2_dbl(Word, Vowel, function(.x, .y) {
    file.path("data-raw", "greek-clusters", "initial-clusters-%s.txt") %>%
      sprintf(.y) %>%
      readr::read_lines() %>%
      match(x = .x)
  })) %>%
  dplyr::mutate(On = (Production * 2) - 1) %>%
  dplyr::mutate(Off = Production * 2) %>%
  dplyr::mutate(Onset = purrr::map2_dbl(Sweep, On, function(.x, .y) {
    FindTagTime(.x, .y)
  })) %>%
  dplyr::mutate(Offset = purrr::map2_dbl(Sweep, Off, function(.x, .y) {
    FindTagTime(.x, .y)
  })) %>%
  dplyr::select(-Production, -On, -Off)










devtools::use_data(
  carstens_filters,
  greek_clusters,
  internal = TRUE,
  overwrite = TRUE
)

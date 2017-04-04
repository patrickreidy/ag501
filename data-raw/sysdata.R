
suppressWarnings(library(tidyverse))


# Internalize the filters that come packaged with the Carstens control server.
carstensFilters <- list.files(path = "data-raw/filters", full.names = TRUE) %>%
  as_tibble() %>%
  rename(Path = value) %>%
  mutate(Name = basename(Path) %>%
           stringr::str_replace(pattern = "\\.txt$", replacement = "")) %>%
  mutate(PassbandLimit = Name %>%
           stringr::str_split(pattern = "_") %>%
           map_chr(4) %>%
           as.numeric()) %>%
  mutate(StopbandLimit = Name %>%
           stringr::str_split(pattern = "_") %>%
           map_chr(5) %>%
           as.numeric()) %>%
  mutate(Attenuation = Name %>%
           stringr::str_split(pattern = "_") %>%
           map_chr(6) %>%
           as.numeric()) %>%
  mutate(SamplingRate = Name %>%
           stringr::str_split(pattern = "_") %>%
           map_chr(7) %>%
           as.numeric()) %>%
  mutate(DifferenceOrder = ifelse(
    stringr::str_detect(Name, "acc"), 2,
    ifelse(
      stringr::str_detect(Name, "vel"), 1,
      0
    )
  )) %>%
  mutate(Window = map(Path, ~ read_lines(.x) %>% as.numeric())) %>%
  select(-Path) %>%
  arrange(DifferenceOrder, PassbandLimit, StopbandLimit, SamplingRate)







devtools::use_data(
  carstensFilters,
  internal = TRUE,
  overwrite = TRUE
)

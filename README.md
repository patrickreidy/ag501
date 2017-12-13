# ag501

An R package for reading and analyzing kinematic data recorded from a Carstens
AG501 electromagnetic articulograph (EMA). Metadata and positional data
(x-, y-, and z-dimensions, and rotation) recorded from EMA sensors are
structured as a nested data frame. Functions are provided for filtering
positional data (e.g., low-pass Butterworth filter) or for calculating derived
time-series (e.g., velocity or acceleration profiles from either forward or
central differencing).

## Installation

To install the current development from GitHub:
```r
devtools::install_github("ag501", username = "patrickreidy")
```

## Examples

#### Metadata for Greek clusters
```r
GreekClusters()
# # A tibble: 4 x 9
#   Sweep  Word Cluster Vowel Anchor  Onset Offset
#   <chr> <chr>   <chr> <chr>  <chr>  <dbl>  <dbl>
# 1  0021 skavi      sk     a      v 15.328 16.315
# 2  0022 skavi      sk     a      v 14.310 15.356
# 3  0021 spala      sp     a      l  6.455  7.382
# 4  0022 spala      sp     a      l  8.004  8.947
# # ... with 2 more variables: POS <chr>, TXT <chr>
```

#### Read sweeps for the data set of Greek clusters
```r
sensor_array <- c("", "", "TB", "TD", "TT", "UL", "LL")

sweeps <- ReadSweep(GreekClusters()$POS, GreekClusters()$TXT, sensors = sensor_array)

greek <- dplyr::bind_cols(
  dplyr::select(GreekClusters(), Sweep, Word, Cluster, Onset, Offset),
  dplyr::select(sweeps, SamplingRate, PositionData)
)
```

#### Apply a low-pass 5th-order Butterworth filter to all channels in PositionData
```r
greek %>%
  dplyr::mutate(LowPassData = Butterworth(PositionData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low"))
```

#### Time-slice the low-pass filtered data to the interval 0.5 seconds before and after the word
```r
greek %>%
  dplyr::mutate(LowPassData = Butterworth(PositionData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low")) %>%
  dplyr::mutate(WordData = TimeSlice(LowPassData, from = Onset-0.5, to = Offset+0.5))
```

#### Estimate velocity and acceleration profiles for each channel in WordData
```r
greek %>%
  dplyr::mutate(LowPassData = Butterworth(PositionData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low")) %>%
  dplyr::mutate(WordData = TimeSlice(LowPassData, from = Onset-0.5, to = Offset+0.5)) %>%
  dplyr::mutate(VelocityData = CentralDifference(WordData, n = 1, order = 1)) %>%
  dplyr::mutate(AccelerationData = CentralDifference(WordData, n = 1, order = 2))
```

#### Compute tongue-tip speed in sagittal (xz) plane, which is relevant for /s/
```r
TT <- function(x, ...) {
  UseMethod("TT", x)
}

TT.data.frame <- function(x, samplingRate, ...) {
  .tt <- 
    dplyr::select(x, Time, TTx, TTz) %>%
    dplyr::mutate(TTx_vel = CentralDifference(TTx, samplingRate = samplingRate)) %>%
    dplyr::mutate(TTz_vel = CentralDifference(TTz, samplingRate = samplingRate)) %>%
    dplyr::mutate(TTxz_spd = sqrt(TTx_vel^2 + TTz_vel^2))
  return(.tt)
}

TT.list <- function(x, samplingRate, ...) {
  .tt <- purrr::map2(x, samplingRate, TT)
  return(.tt)
}

tt<- greek %>%
  dplyr::mutate(LowPassData = Butterworth(PositionData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low")) %>%
  dplyr::mutate(WordData = TimeSlice(LowPassData, from = Onset-0.5, to = Offset+0.5)) %>%
  dplyr::mutate(TT = TT(WordData, SamplingRate))
```

#### Compute lip aperture in vertical dimension, sagittal plane, and three dimensions
``` r
LA <- function(x, ...) {
  UseMethod("LA", x)
}

LA.data.frame <- function(x, ...) {
  .la <-
    dplyr::select(x, Time, ULx, ULy, ULz, LLx, LLy, LLz) %>%
    dplyr::mutate(LAz = sqrt((ULz-LLz)^2),
                  LAxz = sqrt((ULx-LLx)^2 + (ULz-LLz)^2),
                  LAxyz = sqrt((ULx-LLx)^2 + (ULy-LLy)^2 + (ULz-LLz)^2))
  return(.la)
}

LA.list <- function(x, ...) {
  .la <- purrr::map(x, LA)
  return(.la)
}

la <- greek %>%
  dplyr::filter(Cluster == "sp") %>%
  dplyr::mutate(LowPassData = Butterworth(PositionData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low")) %>%
  dplyr::mutate(WordData = TimeSlice(LowPassData, from = Onset-0.5, to = Offset+0.5)) %>%
  dplyr::mutate(LA = LA(WordData))
```


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
# This sensor array utilizes 7 of the 8 sockets available in a Carstens sensin box:
# 1. HL: sensor attached to the left side of the head-correction glasses
# 2. HR: sensor attached to the right side of the head-correction glasses
# 3. TB: sensor attached to the tongue back (body)
# 4. TD: sensor attached to the tongue dorsum
# 5. TT: sensor attached to the tongue tip
# 6. UL: sensor attached to the upper lip
# 7. LL: sensor attached to the lower lip
# Note: the sensors in sockets 3 through 7 (i.e., TB through LL) were all attached
#       within the midsagittal plane, as best as possible.
sensor_array <- c("HL", "HR", "TB", "TD", "TT", "UL", "LL")

sweeps <- ReadSweep(GreekClusters()$POS, GreekClusters()$TXT, sensors = sensor_array)

greek <- dplyr::bind_cols(
  dplyr::select(GreekClusters(), Sweep, Word, Cluster, Onset, Offset),
  dplyr::select(sweeps, SamplingRate, SensorData)
)
```

#### Apply a low-pass 5th-order Butterworth filter to all channels in SensorData
```r
greek %>%
  dplyr::mutate(LowPassData = Butterworth(SensorData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low"))
# # A tibble: 4 x 8
#   Sweep  Word Cluster  Onset Offset SamplingRate            SensorData           LowPassData
#   <chr> <chr>   <chr>  <dbl>  <dbl>        <dbl>                <list>                <list>
# 1  0021 skavi      sk 15.328 16.315          250 <tibble [5,021 x 49]> <tibble [5,021 x 49]>
# 2  0022 skavi      sk 14.310 15.356          250 <tibble [4,740 x 49]> <tibble [4,740 x 49]>
# 3  0021 spala      sp  6.455  7.382          250 <tibble [5,021 x 49]> <tibble [5,021 x 49]>
# 4  0022 spala      sp  8.004  8.947          250 <tibble [4,740 x 49]> <tibble [4,740 x 49]>
```

#### Time-slice the low-pass filtered data to the interval 0.25 seconds before and after the word
```r
greek %>%
  dplyr::mutate(LowPassData = Butterworth(SensorData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low")) %>%
  dplyr::mutate(WordData = TimeSlice(LowPassData, from = Onset-0.25, to = Offset+0.25))
# # A tibble: 4 x 9
#   Sweep  Word ...            SensorData           LowPassData            WordData
#   <chr> <chr> ...                <list>                <list>              <list>
# 1  0021 skavi ... <tibble [5,021 x 49]> <tibble [5,021 x 49]> <tibble [372 x 49]>
# 2  0022 skavi ... <tibble [4,740 x 49]> <tibble [4,740 x 49]> <tibble [387 x 49]>
# 3  0021 spala ... <tibble [5,021 x 49]> <tibble [5,021 x 49]> <tibble [356 x 49]>
# 4  0022 spala ... <tibble [4,740 x 49]> <tibble [4,740 x 49]> <tibble [361 x 49]>
```

#### Estimate velocity and acceleration profiles for each channel in WordData
```r
greek %>%
  dplyr::mutate(LowPassData = Butterworth(SensorData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low")) %>%
  dplyr::mutate(WordData = TimeSlice(LowPassData, from = Onset-0.25, to = Offset+0.25)) %>%
  dplyr::mutate(VelocityData = CentralDifference(WordData, n = 1, order = 1)) %>%
  dplyr::mutate(AccelerationData = CentralDifference(WordData, n = 1, order = 2))
# # A tibble: 4 x 11
#   Sweep  Word ...            WordData        VelocityData    AccelerationData
#   <chr> <chr> ...              <list>              <list>              <list>
# 1  0021 skavi ... <tibble [372 x 49]> <tibble [372 x 49]> <tibble [372 x 49]>
# 2  0022 skavi ... <tibble [387 x 49]> <tibble [387 x 49]> <tibble [387 x 49]>
# 3  0021 spala ... <tibble [356 x 49]> <tibble [356 x 49]> <tibble [356 x 49]>
# 4  0022 spala ... <tibble [361 x 49]> <tibble [361 x 49]> <tibble [361 x 49]>
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

greek %>%
  dplyr::filter(Cluster == "sp") %>%
  dplyr::mutate(LowPassData = Butterworth(SensorData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low")) %>%
  dplyr::mutate(WordData = TimeSlice(LowPassData, from = Onset-0.25, to = Offset+0.25)) %>%
  dplyr::mutate(LA = LA(WordData))
# # A tibble: 2 x 10
#   Sweep  Word ...            WordData                  LA
#   <chr> <chr> ...              <list>              <list>
# 1  0021 spala ... <tibble [356 x 49]> <tibble [356 x 10]>
# 2  0022 spala ... <tibble [361 x 49]> <tibble [361 x 10]>
```

#### Compute tongue-tip speed in sagittal (xz) plane, which is relevant for /s/
```r
TT <- function(x, ...) {
  UseMethod("TT", x)
}

TT.data.frame <- function(x, samplingRate, ...) {
  .tt <- 
    dplyr::select(x, Time, TTx, TTy, TTz) %>%
    dplyr::mutate(TTx_vel = CentralDifference(TTx, samplingRate = samplingRate)) %>%
    dplyr::mutate(TTy_vel = CentralDifference(TTy, samplingRate = samplingRate)) %>%
    dplyr::mutate(TTz_vel = CentralDifference(TTz, samplingRate = samplingRate)) %>%
    dplyr::mutate(TTz_spd = sqrt(TTz_vel^2)) %>%
    dplyr::mutate(TTxz_spd = sqrt(TTx_vel^2 + TTz_vel^2)) %>%
    dplyr::mutate(TTxyz_spd = sqrt(TTx_vel^2 + TTy_vel^2 + TTz_vel^2))
  return(.tt)
}

TT.list <- function(x, samplingRate, ...) {
  .tt <- purrr::map2(x, samplingRate, TT)
  return(.tt)
}

greek %>%
  dplyr::mutate(LowPassData = Butterworth(SensorData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low")) %>%
  dplyr::mutate(WordData = TimeSlice(LowPassData, from = Onset-0.25, to = Offset+0.25)) %>%
  dplyr::mutate(TT = TT(WordData, SamplingRate))
# # A tibble: 4 x 10
#   Sweep  Word ...            WordData                  TT
#   <chr> <chr> ...              <list>              <list>
# 1  0021 skavi ... <tibble [372 x 49]> <tibble [372 x 10]>
# 2  0022 skavi ... <tibble [387 x 49]> <tibble [387 x 10]>
# 3  0021 spala ... <tibble [356 x 49]> <tibble [356 x 10]>
# 4  0022 spala ... <tibble [361 x 49]> <tibble [361 x 10]>
```

#### Find gestural landmarks for movement of the TT sensor
```r
# Landmarks in the sagittal plane only
greek %>%
  dplyr::mutate(LowPassData = Butterworth(SensorData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low")) %>%
  dplyr::mutate(WordData = TimeSlice(LowPassData, from = Onset-0.25, to = Offset+0.25)) %>%
  dplyr::mutate(TT = TT(WordData, SamplingRate)) %>%
  dplyr::mutate(LandmarksTT = FindLandmarks(TT, channels = "TTxz_spd", onsetNear = Onset)) %>%
  dplyr::select(Sweep, Word, Onset, SamplingRate, LandmarksTT) %>%
  tidyr::unnest()
# # A tibble: 16 x 8
#    Sweep  Word  Onset SamplingRate  Channel      Landmark   Time    Value
#    <chr> <chr>  <dbl>        <dbl>    <chr>         <chr>  <dbl>    <dbl>
#  1  0021 skavi 15.328          250 TTxz_spd  GestureOnset 15.328 53.24821
#  2  0021 skavi 15.328          250 TTxz_spd   TargetOnset 15.472 37.50784
#  3  0021 skavi 15.328          250 TTxz_spd  TargetOffset 15.500 26.18081
#  4  0021 skavi 15.328          250 TTxz_spd GestureOffset 15.696 31.36240
#  5  0022 skavi 14.310          250 TTxz_spd  GestureOnset 14.312 45.83411
#  6  0022 skavi 14.310          250 TTxz_spd   TargetOnset 14.476 32.10593
#  7  0022 skavi 14.310          250 TTxz_spd  TargetOffset 14.512 30.39915
#  8  0022 skavi 14.310          250 TTxz_spd GestureOffset 14.720 32.50779
#  9  0021 spala  6.455          250 TTxz_spd  GestureOnset  6.464 54.75046
# 10  0021 spala  6.455          250 TTxz_spd   TargetOnset  6.596 42.25295
# 11  0021 spala  6.455          250 TTxz_spd  TargetOffset  6.624 38.85320
# 12  0021 spala  6.455          250 TTxz_spd GestureOffset  6.788 36.93888
# 13  0022 spala  8.004          250 TTxz_spd  GestureOnset  7.996 62.27941
# 14  0022 spala  8.004          250 TTxz_spd   TargetOnset  8.140 36.79614
# 15  0022 spala  8.004          250 TTxz_spd  TargetOffset  8.236 30.07731
# 16  0022 spala  8.004          250 TTxz_spd GestureOffset  8.412 56.81165

# Landmarks in vertical dimension, sagittal plane, and three-dimensions
greek %>%
  dplyr::mutate(LowPassData = Butterworth(SensorData, SamplingRate,
                                          order = 5, cutoffs = 20, type = "low")) %>%
  dplyr::mutate(WordData = TimeSlice(LowPassData, from = Onset-0.25, to = Offset+0.25)) %>%
  dplyr::mutate(TT = TT(WordData, SamplingRate)) %>%
  dplyr::mutate(LandmarksTT = FindLandmarks(TT, matches = "_spd", onsetNear = Onset)) %>%
  dplyr::select(Sweep, Word, Onset, SamplingRate, LandmarksTT) %>%
  tidyr::unnest()
# # A tibble: 48 x 8
#    Sweep  Word  Onset SamplingRate   Channel      Landmark   Time    Value
#    <chr> <chr>  <dbl>        <dbl>     <chr>         <chr>  <dbl>    <dbl>
#  1  0021 skavi 15.328          250   TTz_spd  GestureOnset 15.316 27.60756
#  2  0021 skavi 15.328          250   TTz_spd   TargetOnset 15.476 23.27646
#  3  0021 skavi 15.328          250   TTz_spd  TargetOffset 15.500 23.80788
#  4  0021 skavi 15.328          250   TTz_spd GestureOffset 15.692 22.94749
#  5  0021 skavi 15.328          250  TTxz_spd  GestureOnset 15.328 53.24821
#  6  0021 skavi 15.328          250  TTxz_spd   TargetOnset 15.472 37.50784
#  7  0021 skavi 15.328          250  TTxz_spd  TargetOffset 15.500 26.18081
#  8  0021 skavi 15.328          250  TTxz_spd GestureOffset 15.696 31.36240
#  9  0021 skavi 15.328          250 TTxyz_spd  GestureOnset 15.328 53.24979
# 10  0021 skavi 15.328          250 TTxyz_spd   TargetOnset 15.472 38.08035
# # ... with 38 more rows
```

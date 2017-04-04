# ag501

Analyze kinematic data recorded from a Carstens AG501 electromagnetic
articulograph (EMA). Metadata and positional data recorded from EMA sensors are 
structured as a nested data frame. Functions provided for filtering data (e.g., 
low-pass Butterworth filter) or for calculating derived time-series from the 
position data (e.g., velocity or acceleration from forward differencing or 
central differencing) are defined with data-first argument lists, making them 
compatible with the `tidyverse`. The forward pipe operator (%>%) is reexported
from `magrittr`.

## Installation

To install the current development from GitHub:
```r
devtools::install_github("ag501", username = "patrickreidy")
```

## Examples

#### Read a single sweep
```r
sweep <- ReadSweep(
  sweep = "0021", 
  path = system.file("extdata", package = "ag501"),
  sensors = c("", "", "TB", "TD", "TT", "UL", "LL")
)
```

#### Apply a low-pass Butterworth filter to all channels of position data
```r
sweep %>%
  dplyr::mutate(
    LowpassData = purrr::map2(
      PositionData, SamplingRate,
      Butterworth,
      order = 5, cutoffs = 20, type = "low"
    )
  ) %>%
  dplyr::select(
    Sweep, SamplingRate, TimeData, 
    PositionData, LowpassData
  )
```

#### Estimate velocity and acceleration of all channels of position data
```r
sweep %>%
  dplyr::mutate(
    VelocityData = purrr::map2(
      PositionData, SamplingRate,
      CentralDifference,
      n = 1, order = 1, suffix = "_vel"
    )
  ) %>%
  dplyr::mutate(
    AccelerationData = purrr::map2(
      PositionData, SamplingRate,
      CentralDifference,
      n = 1, order = 2, suffix = "_acc"
    )
  ) %>%
  dplyr::select(
    Sweep, SamplingRate, TimeData,
    PositionData, VelocityData, AccelerationData
  )
```

#### Estimate speed of tongue-tip sensor (TT) in the sagittal (xz) plane
```r
sweep %>%
  dplyr::mutate(
    TongueTipData = purrr::pmap(
      list(TimeData, PositionData, SamplingRate),
      function(td, pd, sr) {
        dplyr::bind_cols(td, pd) %>%
          dplyr::select(Time, TTx, TTz) %>%
          dplyr::mutate(TTx_vel = CentralDifference(TTx, sr)) %>%
          dplyr::mutate(TTz_vel = CentralDifference(TTz, sr)) %>%
          dplyr::mutate(TTxz_spd = sqrt(TTx_vel^2 + TTz_vel^2))
     }
    )
  ) %>%
  dplyr::select(
    Sweep, SamplingRate, TimeData,
    PositionData, TongueTipData
  )
```

#### Low-pass filter position data; compute lip aperature time-series
```r
sweep %>%
  dplyr::mutate(
    LowpassData = purrr::map2(
      PositionData, SamplingRate,
      Butterworth,
      order = 5, cutoffs = 20, type = "low"
    )
  ) %>%
  dplyr::mutate(
    LipApertureData = purrr::pmap(
      list(TimeData, LowpassData, SamplingRate),
      function(td, pd, sr) {
        dplyr::bind_cols(td, pd) %>%
        dplyr::select(Time, ULx, ULz, LLx, LLz) %>%
        dplyr::mutate(APxz = sqrt((ULx-LLx)^2 + (ULz-LLz)^2)) %>%
        dplyr::mutate(APxz_vel = CentralDifference(APxz, sr))
      }
    )
  ) %>%
  dplyr::select(
    Sweep, SamplingRate, TimeData, 
    PositionData, LowpassData, LipApertureData
  )
```

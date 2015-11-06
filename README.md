## This is a R package to help working with biologging data, focused on diving predators.

Most of the functions in this package were coded to work with Southern elephant seal data but they must be useful for other diving predator such as penguins, turtles, other seals etc... Let me know your study species and the settings you use.

### Motivation

It seems to me that code sharing is a key to **reproductible research** and **collaboration among researchers**. Distributing analysis tools is also a way to **bring biologging data to new people** as laboratories are increasingly sharing their data.

This open source package aims at sharing my work as well as at creating a place where other people working on diving predators can help me to improve it, contribute with their own work and maybe interact with enthusiastic colleagues from all over the world.

--------

## The package includes documented functions for many tasks

### Process raw data

With functions to correct common issues with TDR (Time-Depth Recorder) data and to process acceleration/magnetometry data.

  * Detect and fix gaps in **Time-Depth profiles**
  * Fix drift in depth records over time.
  * Detect **Prey Catch Attempts**

### Identify specific parts in a dataset 

  * Dive identification
  * **Bottom of dives** delimitation (4 different methods are implemented)
  * **Dift dives** detection

### Brokenstick Models (BSM)

These tools allow to quickly build and manipulate BSM thus enabling working with 
old datasets (6 points BSM dive profiles) and newer ones (high sampling frequency).

  * **Fit brokenstick models** with a fixed number of points or according to a cost function
  * Plot BSM objects
  * Extract BSM parameters (**breakpoints**, **slopes**, intercept)
  * Extract fitted values or residuals, make predictions
  * Compute **dive zone index**

### Facilitate data manipulation, extraction and visualization of animal behavior 

With methods for plot, brokenstick models, tdr tables, dive statistics table and 
functionnal programming utilities.
With diving behavior indices such as **number of wiggles**, **time-at-depth**, 
**sinuosity**

  * Dives can be subseted using their ID.
  * Symbols can be used to extract specific parts of a dive cycle represented as `"~!_/-"`. `"!_/"` select dives, `"!"` select the descents, `"_"` select the bottoms, `"/"` the select ascents, `"!/"` select descents and asents nd merge them into a single table while `"!&/"` select descents and asents but keep them separated in distinct tables. Brakets and OR operators are also implemented. Surface periods are designated by `"~"` and `"-"` depending on their position with reference to dives.

**Few examples** (with the dataset included in the package) in order show how you may use the functions provided by the `rbl` package:

Basic usage of `tdrply`: 

`tdrply(f = function to apply, cl = TDR columns, ty = dive cycle part, no = dive number, ... = other arguments to f)`

```S
data(exses) # Load the package example dataset
ind(ses)    # Set it as default individual to avoid writing its name everytime

# Extract bottom of the dive number 65 of exses
dv_btt <- tdrply(identity, ty = "_", no = 65)

# Maximum depth of the dive number 400
tdrply(max, "depth", ty = "!_/", no = 400, na.rm = TRUE)
# !_/#400 
# 428.9812 
```

Usage of `ty` argument
```S
# Average pitch angles of (descent + ascent) phases of dives no 400 to 405
# ty = "!/" returns a single vectors with results computed on (descent + ascent) 
tdrply(mean, "pitch", ty = "!/", no = 400:405)
#      !/#400      !/#401      !/#402      !/#403      !/#404      !/#405 
# -0.29398292 -0.09295794 -0.12892324 -0.23085102  0.02746734 -0.07436099 

# Average pitch angles of descent and ascent phases of dives no 400 to 405
# ty = "!&/" returns a single vectors with results computed on (descent) + (ascent) 
tdrply(mean, "pitch", ty = "!&/", no = 400:405)
#     !#400      /#400      !#401      /#401      !#402      /#402 
# -1.0048280  0.5184114 -1.1039254  0.8592323 -1.3239197  0.6006535 
#      !#403      /#403      !#404      /#404      !#405      /#405 
# -0.9891836  0.5491482 -1.0761095  0.5225652 -0.7189206  0.6852277 

# Average pitch angles of descent and ascent phases of dives no 400 to 405
# ty = c("!", "/") returns a list of two vectors, results for 1-descents and 2-ascents
tdrply(mean, "pitch", ty = c("!", "/"), no = 400:405)
# $dsc
#      !#400      !#401      !#402      !#403      !#404      !#405 
# -1.0048280 -1.1039254 -1.3239197 -0.9891836 -1.0761095 -0.7189206 
# 
# $asc
#     /#400     /#401     /#402     /#403     /#404     /#405 
# 0.5184114 0.8592323 0.6006535 0.5491482 0.5225652 0.6852277 
```

Quick fitting and plotting of brokenstick model
```S
# 6 breakpoints BSM (and BSM plot) of dive number 477 to 480
bsm <- tdrply(brokenstick, c("time", "depth"), ty = "!_/", no = 477:480, npts = 6)
par(mfrow = c(2, 2), mar = c(4, 4, .2, .2))
lapply(bsm, plot, 
      data = TRUE, # Add the high sampling frequency dive profile
      enumerate = TRUE # Enumerate the breakpoints
)
```

![Brokenstick models (dives 477 to 480)](http://oi67.tinypic.com/4rdrbk.jpg "Brokenstick models (dives 477 to 480)")

Count and visualize wiggles in a dive
```S
# Wiggles in dive no 480
par(mfrow = c(1, 1))
tdrply(wiggles, c("time", "depth"), ty = "!_/", no = 480, 
      plt = TRUE # to make a plot of the identified wiggles
      )
# !_/#480 
#      18 
```

![Wiggles in dive 480](http://oi64.tinypic.com/2ic5ci1.jpg "Wiggles in dive 480")

Plot dive profiles and compute behavioral variables
```S
# Profile with Prey Catch Attempts of dive number 480
# The third variable "is_pca" is color coded TRUE is indicated by red color.
tdrply(plot, c("time", "depth", "is_pca"), ty = "!_/", no = 480, 
       main = "plot.tdr example", pch = 20 # Graphical parameters
       )

# Count the number of Prey Catch Attempts
# A continuous succession of TRUE is interpreted as a single PCA
tdrply(pca_count, "is_pca", ty = "_", no = 480)
# _#480 
#    10 
   
# TAD of dive number 480
tdrply(time_at_depth, c("time", "depth"), ty = "!_/", no = 480)
#   !_/#480 
# 0.7645625

# Sinuosity of the bottom of the dive number 480
tdrply(sinuosity, c("time", "depth"), ty = "_", no = 480)
#    _#480 
# 6.452312 
```

![Dive profile with PCA](http://oi66.tinypic.com/2nkkmb.jpg "Dive profile with PCA")

--------

## Package binaries are available at 

https://github.com/SESman/rbl/releases.

```S
#### Install using dowloaded binaries ####
install.packages("path/to/the/downloaded/binary.extension", repos = NULL)

#### Install with devtools ####
install.packages("devtools")
require("devtools")
install_github("rbl", "SESman")

#### Overview of package ####
require(rbl) # Make sure you have installed the dependencies (e.g. RcppRoll)
help(package = "rbl")
```

## Please feel free to contribute to the project by adding your code, reporting the bugs or acknowledging the work.

#' Reciprocal function of base::log10
#' 
#' @param x numeric vector
#' @export
#' @keywords internal light
exp10 <- function(x) exp(x * log(10))

#' Convert Wildlife Computers light values from/to linear W/cm^2 units
#' 
#' Use the equation provided by the manufacturer, Wildlife Computers.
#' 
#' @param x raw/transformed sensor readings
#' @export
#' @keywords light
SI_light <- function(x) {
  10^((x - 250) / 20)
}

#' @rdname SI_light
#' @export
#' @keywords light
WC_light <- function(x) {
  20 * log10(x) + 250
}

#' BioPIC QR decomposition on a TDR sample
#' @param x A data subset of fixed width of 11 seconds/lines.
#' @param lightSI.nm The name of the TDR column with light values in W/cm^2.
#' @references 
#' Vacquie-Garcia, J., Royer, F., Dragon, A.-C., Viviant, M., Bailleul, F. 
#' & Guinet, C. (2012) Foraging in the Darkness of the Southern Ocean: 
#' Influence of Bioluminescence on a Deep Diving Predator. PLoS ONE, 7, e43565.
#' @keywords internal
#' @export
biopic.qr <- function(x, lightSI.nm = "light_si") {
  # A = QR. We use qr.solve() to find R given A and Q.
  Amat <- log10(x[ , lightSI.nm])
  # Build the Q matrix (w rows X length(R) columns)
  Qmat <- matrix(c(
    0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,  # log(It) = log(Ia)
    1, 0, 1,                                      # log(It) = log(alpha) + log(Ia)
    1, 1, 1, 1, 2, 1, 1, 3, 1, 1, 4, 1, 1, 5, 1), # log(It) = log(alpha) -K*t + log(Ia)
    nrow = 11, byrow = TRUE)
  Rmat <- try(qr.solve(Qmat, Amat), silent = TRUE)
  if (is.error(Rmat)) {
    Rmat <-rep(NA, 3)
  } else {
    Rmat[c(1, 3)] <- exp10(Rmat[c(1, 3)])
  }
  Rmat # alpha, lessK & Ia
}

#' BioPIC: bioluminescent event detection tool
#' 
#' Adaptation of BioPIC method with 1 Hz sampling frequency datasets and 
#' of post-peak attenuation (related to sensor and not depth)
#' 
#' @param x a TDR data sample
#' @param nms The names of the output variables.
#' @inheritParams biopic.qr
#' @export
#' @keywords light
#' @references 
#' Vacquie-Garcia, J., Royer, F., Dragon, A.-C., Viviant, M., Bailleul, F. 
#' & Guinet, C. (2012) Foraging in the Darkness of the Southern Ocean: 
#' Influence of Bioluminescence on a Deep Diving Predator. PLoS ONE, 7, e43565.
#' @examples 
#' data(exses)
#' dv <- tdrply(identity, no = 100, obj = exses)[[1]]
#' dv$light_si <- SI_light(dv$light)
#' biolum_tbl <- BioPIC(dv)
BioPIC <- function(x, lightSI.nm = "light_si", nms = c("alpha", "lessK", "Ia")) {
  x[ , nms] <-  NA
  nr <- nrow(x)
  for (ii in seq(6, nr - 5)) {
    x[ii, nms] <- biopic.qr(x[seq(ii - 5, ii + 5), ], lightSI.nm)
  }
  attr(x, "biopic") <- list("lightSI.nm" = lightSI.nm, "bioPIC.nms" = nms)
  x
}

#' Identify potential bioluminescence emission events
#' 
#' @param x a data.frame such as returned by \code{\link{BioPIC}} output
#' @param sensitivity A threshold to decide if a signal peak is high enought. See 
#' details.
#' @param max_light Max light level (W/cm^2) where a event can be detected. Set to 
#' NULL to disable.
#' @details Sensor sensitivity threshold: The minimum ratio between two measures 
#' to be considered significantly different. In laboratory experiments, the 
#' sensors measured light intensity +- 2 units while submitted to a constant 
#' light intensity. Translating the sensor log-scale units to SI linear scale 
#' units this accuracy measure translates into a minimum ratio of 1.26
#' @export
#' @keywords light
#' @examples 
#' \dontrun{
#' data(exses)
#' exses$tdr$light_si <- SI_light(exses$tdr$light)
#' biolum_tbl <- tdrply(BioPIC, obj = exses)
#' ble_tbl <- Reduce(rbind, lapply(biolum_tbl, biolum_events))
#' }
biolum_events <- function(x, sensitivity = 1.26, max_light = exp(-20)) {
  # Retrieve info
  alpha <- attr(x, "biopic")$bioPIC.nms[1]
  light <- attr(x, "biopic")$lightSI.nm[1]
  
  # Increasing light periods
  is_potble <- is.finite(x[ , alpha]) & x[ , alpha] > 1 # Alpha is valid and > 1
  potble <- per(is_potble)
  potble <- potble[potble$value %in% TRUE, ]
  
  # Compute peak light ratio
  potble$st_idx <- potble$st_idx
  potble$start_time  <- x[potble$st_idx, 1]
  potble$end_time    <- x[potble$ed_idx, 1]
  potble$start_light <- x[potble$st_idx, light]
  peakmax_idx <- mapply(function(st, ed) which.max(x[st:ed, light]), potble$st_idx, potble$ed_idx)
  peakmax_idx <- peakmax_idx + potble$st_idx - 1
  potble$peakmax_time  <- x[peakmax_idx, 1]
  potble$peakmax_light <- x[peakmax_idx, light]
  potble$peak_ratio <- potble$peakmax_light / potble$start_light
  
  # Filter
  cnd <- potble$peak_ratio >= sensitivity
  cnd <- "if"(is.null(max_light), cnd, cnd  & potble$peakmax_light <= max_light)
  potble[cnd, -(1:4)]
}

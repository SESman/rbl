#' Southern ocean climatologic fronts
#' 
#' This dataset contains the coordinates of the the major front of the 
#' southern ocean: SACCF (Southern Antarctic Circumpolar Current Front), 
#' PF (Polar Front), SAF (Sub-Antarctic Front), 
#' SSTF (Southern Sub Tropical Front) and NSTF.
#' 
#' The SSTF position is defined by the intersection
#' between the 11 deg. C isotherm and the 150 m isobath. It represents the limits between the northern
#' warm and salty waters (Atlantic, Indian and Pacific oceans) and the colder and less salty
#' water from the sub-antarctic area. It is considered as the northern limit of the antarctic circumpolar
#' current (but not of the Southern Ocean which has as been arbitrarily defined to the 40 deg. S
#' 
#' The SAF position is defined by the maximum meridian temperature
#' gradient between 3 deg. and 8 deg. C at 300 m. It is associated with strong currents. Around
#' Kerguelen it deviates to the north, thus reducing the extent of the sub-antarctic area and increasing
#' that of the polar-frontal area.
#' 
#' The PF position is defined as the northern limit of the minimum subsurface
#' temperature lower than 2 deg. C. It represents the lower trace of the winter-mixed layer, which
#' the upper part warms up during summer. It is associated with strong currents.
#' 
#' The SACCF front is defined as the southern limit of the Antarctic circumpolar 
#' current. It represents the border between the antarctic
#' (to the north) and continental (to the south) areas.
#' 
#' Content:
#' \itemize{
#'   \item Lon Longitude
#'   \item Lat Latitude
#'   \item name The front name.
#' }
#' 
#' @name sofronts
#' @usage data(sofronts)
#' @format A data frame of 9212 rows and 3 variables.
#' @keywords datasets
#' @docType data
#' @seealso \code{\link{front}}
#' @references
#' Belkin, I.M. (1988) Main hydrological features of the Central South Pacific, in:  Ecosystems  of the 
#' Subantarctic Zone of the Pacific Ocean, edited by M.E. Vinogradov and M.V. Flint, Nauka, 
#' Moscow, 21 28 [Translated as "Pacific Subantarctic Ecosystems", pp.12-17, New Zealand Translation Centre Ltd., Wellington, 1996].
#' 
#' Belkin, I.M. (1993) Frontal structure of the South Atlantic, in: Pelagic Ecosystems of the 
#' Southern Ocean, edited by N.M. Voronina, pp. 40 53 (in Russian), Nauka, Moscow.
#' 
#' Belkin, I.M., and A.L. Gordon (1996) Southern Ocean fronts from the Greenwich 
#' meridian to Tasmania, J. Geophys. Res., 101(C2), 3675-3696.
#' 
#' Wessel, P., and W. H. F. Smith (1996) A Global Self-consistent, Hierarchical, 
#' High-resolution Shoreline Database, J. Geophys. Res., 101(B4), 8741-8743. 
NULL

#' Bathymetry of the Kerguelen area
#' 
#' This dataset is a imported shapefile containing the Hi-Res. isobath of the 
#' Kerguelen area. Available isobaths from -200 m to -4400 meters by 200 m steps.
#' 
#' @name kerbathy
#' @usage data(kerbathy)
#' @format A SpatialLinesDataFrame (S4 object of \code{sp} package)
#' @keywords datasets
#' @docType data
#' @seealso \code{\link{isobath}}
#' @references
#' \url{https://www.ga.gov.au/products/servlet/controller?event=GEOCAT_DETAILS&catno=71552}
NULL

#' An tiny dataset (about 10 days) of an SES equipped with an accelerometer
#' 
#' @name exses
#' @usage data(exses)
#' @format An object of class \code{'ses'}.
#' @keywords datasets
#' @docType data
#' @seealso \code{\link{as.ses}}
NULL

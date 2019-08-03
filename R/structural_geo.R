#' Converts from ENU (east, north, up) to NED (north, east, down) coordinates
#' @param enu_vec a vector in ENU coordinates
#' @return a vector in NED coordinates
#' @export
enu2ned <- function(enu_vec){
  ned_vec <- c(enu_vec[2], enu_vec[1], -enu_vec[3])
  names(ned_vec) <- c('north','east','down')
  return (ned_vec)
}

#' Converts from NED (north, east, down) to ENU (east, north, up) coordinates
#' @param ned_vec a vector or matrix in NED coordinates, with NED in the columns
#' @return a vector in ENU coordinates
#' @export
ned2enu <- function(ned_vec){
  enu_vec <- c(ned_vec[2], ned_vec[1], -ned_vec[3])
  names(enu_vec) <- c('east','north','up')
  return (enu_vec)
}

#' Converts from cartesian to spherical coordinates
#' @param ned_vec this is a vector of direction cosines (i.e. unit vector) in NED coordinates
#' @return an array of the strike and dip (or trend and plunge) in the lower hemispher projection, in radians
#' @export
cart2sph <- function(ned_vec){
  dip <- asin(ned_vec[3])
  if (ned_vec[1] == 0){
    if (ned_vec[2] < 0){
      strike <- 3/2*pi # west strike
    } else {
      strike <- pi/2 # east strike
    }
  } else {
    strike <- atan(ned_vec[2]/ned_vec[1])
    if (ned_vec[1] < 0){strike <- strike + pi}
  }
  # constrain azimuth between 0 and 2pi
  strike <- strike - floor(strike / (2*pi)) * 2*pi

  sd_vec <- c(strike,dip)
  names(sd_vec) <- c('strike','dip')

  return (sd_vec)
}

#' Converts from spherical to cartesian coordinates
#' @param strike_or_trend this is a scalar of the trend or strike in radians
#' @param dip_or_plunge this is a scalar of the plunge or dip in radians
#' @param plane this is a boolean that selects whether the strike/dip of a plane (TRUE) or trend/plunge of a line (FALSE) are analyzed
#' @return an array of the north,east,down (NED) direction cosines (i.e. unit vector)
#' @export
sph2cart <- function(strike_or_trend = pi,dip_or_plunge = pi/2, plane = TRUE){
  if (plane == TRUE){
    ned_vec <- c(sin(dip_or_plunge) * sin(strike_or_trend),
                 -sin(dip_or_plunge) * cos(strike_or_trend),
                 cos(dip_or_plunge))
  } else {
    ned_vec <- c(cos(dip_or_plunge) * cos(strike_or_trend),
                 cos(dip_or_plunge) * sin(strike_or_trend),
                 sin(dip_or_plunge))
  }

  names(ned_vec) <- c('north','east','down')
  return (ned_vec)
}

#' Converts from the trend/plunge of a pole to the strike/dip of a plane
#' @param trend this is a scalar of the pole trend in radians
#' @param plunge this is a scalar of the pole plunge in radians
#' @return an array of the strike and dip of the plane in radians
#' @export
pole2plane <- function(trend = pi,plunge = pi/2){
  if (plunge >= 0){
    dip <- pi/2 - plunge
    strike <- trend + pi/2
  } else {
    dip <- pi/2 + plunge
    strike <- trend - pi/2
  }

  # constrain azimuth between 0 and 2pi
  strike <- strike - floor(strike / (2*pi)) * 2*pi

  sd_vec <- c(strike,dip)
  names(sd_vec) <- c('strike','dip')

  return (sd_vec)
}

#' Converts from the strike/dip of a plane to the trend/plunge of a pole
#' @param strike this is a scalar of the plane strike in radians
#' @param dip this is a scalar of the plane dip in radians
#' @return an array of the strike and dip of the plane in radians
#' @export
plane2pole <- function(strike = pi,dip = pi/2){
  u_vec <- sph2cart(strike,dip,plane = TRUE)
  sd_vec <- cart2sph(u_vec)

  return (sd_vec)
}

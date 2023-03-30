#' Function to compute derivative of likelihood with respect to the dispersion parameter \eqn{\phi}
#' 
#' Internal Function
#' 
#' @param y observed response
#' @param xi index parameter
#' @param mu mean vector
#' @param phi dispersion vector 
#' @keywords 
#' @import stats tweedie 
#' @export
#' @examples 

delphi <- function(y = NULL,
                   xi = NULL,
                   mu = NULL,
                   phi = NULL,
                   h = 1e-10){
  dtw.h = dtweedie(y = y, xi = xi, mu = mu, phi = phi + h)
  dtw.h[dtw.h == 0] = 1e-300
  dtw = dtweedie(y = y, xi = xi, mu = mu, phi = phi)
  dtw[dtw == 0] = 1e-300
  d1phi = (log(dtw.h) - log(dtw))/h
  
  return(d1phi)
}

#' Deviance function for Tweedie Compound Poisson-gamma distribution
#'
#' Internal Function
#'
#' @param y observed response
#' @param wts optional weights (or scales) for y
#' @param mu mean vector
#' @param phi dispersion vector
#' @param xi index parameter
#' @keywords
#' @import stats
#' @importFrom
#' @export
#' @examples
#'
dev.tw <- function(y = NULL,
                   wts = NULL,
                   mu = NULL,
                   phi = NULL,
                   xi = NULL){
  if(is.null(wts)){
    (y * mu^(1 - xi)/(1 - xi) - mu^(2 - xi)/(2 - xi))/phi
  }else{
    2 * wts * (y^(2 - xi) / ((1 - xi) * (2 - xi)) - y * mu^(1 - xi)/(1 - xi) + mu^(2 - xi)/(2 - xi))/phi
  }
}

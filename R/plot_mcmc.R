#' MCMC Diagnostic Plots for samplers
#'
#' Internal Function
#'
#' @param samples MCMC samples
#' @param true true values, optional
#' @param cnames parameter names
#' @keywords
#' @import stats
#' @export
#' @examples

plot_mcmc <- function(samples = NULL,
                      true = NULL,
                      col = "darkgreen",
                      acf.lag.max = 100,
                      cnames = NULL){
  nc = ncol(samples)
  if(is.null(nc)) nc = 1
  N = nrow(samples)
  if(is.null(cnames)) cnames = colnames(samples)

  if(nc == 1){
    ts.plot(samples, main = paste("Trace of", cnames, sep = " "), ylab = "", xlab = "Iterations", col="blue")
    abline(h = median(samples), lty="dotted", col = col)
    if(!is.null(true)) abline(h = true, lty="dotted", col = "darkred")

    plot(density(samples), main = paste("Density of", cnames, sep = " "), ylab = "",
         xlab = paste("N =", N, "Bandwidth =", sprintf("%.3f", round(density(samples)$bw,3)), sep = " "), col="blue")
    rug(jitter(samples))
    abline(v = median(samples), lty="dotted", col = col)
    if(!is.null(true)) abline(v = true, lty="dotted", col = "darkred")

    acf(samples, lag.max = acf.lag.max, main = paste("ACF for", cnames, sep = " "))
  }else{
    if(sum(par()$mfrow==c(4, 3)) < 2) par(mfrow=c(4, 3))
    for(i in 1:nc){
      ts.plot(samples[,i], main = paste("Trace of", cnames[i], sep = " "), ylab = "", xlab = "Iterations", col="blue")
      abline(h = median(samples[,i]), lty="dotted", col = col)
      if(!is.null(true)) abline(h = true[i], lty="dotted", col = "darkred")

      plot(density(samples[,i]), main = paste("Density of", cnames[i], sep = " "), ylab = "",
           xlab = paste("N =", N, "Bandwidth =", sprintf("%.3f", round(density(samples[,i])$bw,3)), sep = " "), col="blue")
      rug(jitter(samples[,i]))
      abline(v = median(samples[,i]), lty="dotted", col = col)
      if(!is.null(true)) abline(v = true[i], lty="dotted", col = "darkred")

      acf(samples[,i], lag.max = acf.lag.max, main = paste("ACF for", cnames[i], sep = " "))
    }
  }
}

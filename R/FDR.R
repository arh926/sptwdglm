#' Computes False discovery rate using Bayesian q-values
#'
#' @param chain MCMC chain for beta or gamma in DGLM
#' @param thres threshold value
#' @param FDRc percent at which FDR is computed
#' @keywords
#' @import
#' @importFrom
#' @export
#' @examples
#' \dontrun{
#' rbeta.ssdglm = FDR(mc$beta.mcmc)
#' rgamma.ssdglm = FDR(mc$gamma.mcmc)
#' }
FDR <- function(chain = NULL,
                thres = 0.05,
                FDRc = 0.05){
  mIter = nrow(chain)
  nc = ncol(chain)
  x = sapply(1:nc, function(k) mean(abs(chain[,k]) < thres))
  x2 = cumsum(x[order(x)])/1:nc
  r = rep(NA, nc)
  x.ind = order(x)[which(as.numeric(x2<FDRc)==1)]
  r[x.ind] = x[x.ind]
  r
}

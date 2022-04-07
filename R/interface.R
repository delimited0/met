#' @export
rtmvnorm = function(n, mu, Sigma, lb, ub, A = NULL) {
  
  if (is.null(A)) {
    samples = t(mvrandn(lb, ub, Sigma, n, mu))
  }
  
  return(samples)
}
#' @export
rtmvn = function(n, mu, Sigma, lb, ub, A = NULL) {
  
  if (is.null(A)) {
    samples = t(mvrandn(lb, ub, Sigma, n, mu))
  }
  
  return(samples)
}

#' @param n number of 
#' @export
pmvn = function(mu, Sigma, lb, ub, A = NULL,
                   n = 10000, n_est = 12, type = "qmc") {
  
  # axis aligned problem
  if (is.null(A)) {
    
    lb_center = lb - mu
    ub_center = ub - mu
    
    if (type == "qmc") {
      result = mvNqmc(lb_center, ub_center, Sigma, n, n_est)
    }
    
  }
  
  prob = result$prob
  attr(prob, "error") = result$error
  attr(prob, "relErr") = result$relErr
  attr(prob, "upbnd") = result$upbnd
  
  return(prob)
}
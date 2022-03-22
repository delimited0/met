#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

arma::vec lnNpr(arma::vec a, arma::vec b) {

  int d = a.n_elem;  
  double prob_a;
  double prob_b;
  
  arma::vec probs(d);
  
  for (int i = 0; i < d; i++) {
    
    // case b > a > 0
    if (a(i) > 0) {
      prob_a = R::pnorm(a[i], 0.0, 1.0, 0, 1);
      prob_b = R::pnorm(b[i], 0.0, 1.0, 0, 1);
      probs(i) = prob_a + std::log1p(-std::exp(prob_b - prob_a));
    }
    // case a < b < 0
    else if (b(i) < 0) {
      prob_a = R::pnorm(a[i], 0.0, 1.0, 1, 1);
      prob_b = R::pnorm(b[i], 0.0, 1.0, 1, 1);
      probs(i) = prob_b + log1p(-exp(prob_a - prob_b));
    }
    // case a < 0 < b
    else {
      prob_a = R::pnorm(a[i], 0.0, 1.0, 1, 0);
      prob_b = R::pnorm(b[i], 0.0, 1.0, 0, 0);
      probs[i] = log1p(-prob_a - prob_b);
    }
  }
  
  return probs;
}

double lnNpr(double a, double b) {
  
  double prob_a;
  double prob_b;
  
  double logprob;
  
  // case b > a > 0
  if (a > 0) {
    prob_a = R::pnorm(a, 0.0, 1.0, 0, 1);
    prob_b = R::pnorm(b, 0.0, 1.0, 0, 1);
    logprob = prob_a + std::log1p(-std::exp(prob_b - prob_a));
  }
  // case a < b < 0
  else if (b < 0) {
    prob_a = R::pnorm(a, 0.0, 1.0, 1, 1);
    prob_b = R::pnorm(b, 0.0, 1.0, 1, 1);
    logprob = prob_b + log1p(-exp(prob_a - prob_b));
  }
  // case a < 0 < b
  else {
    prob_a = R::pnorm(a, 0.0, 1.0, 1, 0);
    prob_b = R::pnorm(b, 0.0, 1.0, 0, 0);
    logprob = log1p(-prob_a - prob_b);
  }
  
  return(logprob);
}
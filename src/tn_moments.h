#ifndef TN_MOMENTS_H
#define TN_MOMENTS_H

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

arma::vec lnNpr(arma::vec a, arma::vec b);
double lnNpr(double a, double b);

#endif
#include "RcppArmadillo.h"
#include <boost/math/distributions/normal.hpp>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include "tn_moments.h"

const double INV_THRESHOLD = 0.4;
const double TOL = 2.05;

// accept-reject standard normal truncated to [l, u], l > 0
// Rayleigh proposal distribution
// [[Rcpp::export]]
arma::vec ntail(arma::vec l, arma::vec u) {
  
  int n = l.n_elem;
  
  arma::vec c = arma::square(l) / 2.;
  arma::vec q = arma::expm1(c - arma::pow(u, 2.) / 2.);
  
  arma::vec x(n);
  double prop;
  
  for (int i = 0; i < n; i++) {
    do {
      prop = c(i) - std::log(1 + arma::randu<double>() * q(i));
    } while (std::pow(arma::randu<double>(), 2.) * prop > c(i));
    x(i) = prop;
  }
  
  return arma::sqrt(2 * x);
}

double ntail(const double l, const double u) {
  double c = std::pow(l, 2.) / 2.;
  double q = std::expm1(c - std::pow(u, 2.) / 2.);
  double x;
  do {
    x = c - std::log(1 + arma::randu<double>() * q);
  } while (std::pow(arma::randu<double>(), 2.) * x > c);
  
  return std::sqrt(2. * x);
}

// [[Rcpp::export]]
arma::vec trnd(arma::vec l, arma::vec u) {
  
  int n = l.n_elem;
  arma::vec x(n);
  double prop;
  
  for (int i = 0; i < n; i++) {
    do {
      prop = arma::randn<double>();
    } while (prop < l(i) || prop > u(i));
    x(i) = prop;
  }
  
  return x;
}

double trnd(const double l, const double u) {
  double x;
  do {
    x = arma::randn<double>();
  } while (x < l || x > u);
  
  return x;
}

// [[Rcpp::export]]
double trninv(const double l, const double u) {
  boost::math::normal dist(0.0, 1.0);
  double pl = .5 * std::erfc(-l / std::sqrt(2.));
  double pu = .5 * std::erfc(-u / std::sqrt(2.));
  double x = quantile(dist, pl + (pu - pl) * arma::randu<double>());
  
  return(x);
}

// [[Rcpp::export]]
arma::vec trandn(const arma::vec l, const arma::vec u) {
  int n = l.n_elem;
  arma::vec x(n);
  
  for (int i = 0; i < n; i++) {
    
    // case 1: a < l < u
    if (l(i) > INV_THRESHOLD) {
      x(i) = ntail(l(i), u(i));
    }
    
    // case 2: l < u < -a
    else if (u(i) < -INV_THRESHOLD) {
      x(i) = -ntail(-u(i), -l(i));
    }
    
    else {
      // case 3: abs(u - l) > tol
      if (std::abs(u(i) - l(i)) > TOL) {
        x(i) = trnd(l(i), u(i));
      }
      // case 4: abs(u - l) < tol
      else {
        x(i) = trninv(l(i), u(i));
      }
    }
  }
  
  return(x);
}

double trandn(const double l, const double u) {
  double x;
  
  // case 1: a < l < u
  if (l > INV_THRESHOLD) {
    x = ntail(l, u);
  }
  
  // case 2: l < u < -a
  else if (u < -INV_THRESHOLD) {
    x = -ntail(-u, -l);
  }
  
  else {
    // case 3: abs(u - l) > tol
    if (std::abs(u - l) > TOL) {
      x = trnd(l, u);
    }
    // case 4: abs(u - l) < tol
    else {
      x = trninv(l, u);
    }
  }
  
  return x;
}

/* 
 * n: number of samples
 * mu: optimal tilting parameter
 * psistar: acceptance threshold ?
 * L: permuted covariance cholesky factor
 * l: permuted lower bound
 * u: permuted upper bound
*/
// [[Rcpp::export]]
arma::mat mvn_tilted_accept_reject(int n, arma::vec mu, double psistar, 
                                 arma::mat L, arma::vec l, arma::vec u) {
  
  int d = l.n_elem;
  
  int n_accepted = 0;
  int n_sim = n;
  int n_tot_sim = 0;
  
  arma::mat samples(d, n);
  
  while (n_accepted < n) {
    // Rcpp::Rcout << "accepted: " << n_accepted << std::endl;  
    
    // Rcpp::Rcout << "filling" << std::endl;
    // Rcpp::Rcout << "mu length: " << mu.n_elem << std::endl;
    // Rcpp::Rcout << "d: " << d << std::endl;
    // mu(d-1) = 0; // why do we do this?
    double logpr = 0;  
    
    // Rcpp::Rcout << "Init Z" << std::endl;
    arma::vec Z(d, arma::fill::zeros);
    
    for (int k = 0; k < d; k++) {
      
      // Rcpp::Rcout << "dot product" << std::endl;
      double col = arma::dot(L.submat(k, 0, k, k), Z.rows(0, k));
      
      double tl = l(k) - mu(k) - col;
      double tu = u(k) - mu(k) - col;
      
      Z(k) = mu(k) + trandn(tl, tu);
      
      // Rcpp::Rcout << "computing log prob" << std::endl;
      logpr += lnNpr(tl, tu) + .5 * mu(k) * (mu(k) - mu(k) * Z(k));
    }
    
    // accept reject
    bool accept = arma::randg(arma::distr_param(1., 1.)) > psistar - logpr;
    if (accept) {
      samples.col(n_accepted) = Z;
      // iter++;
      n_accepted++;
    }
  }
  
  return(samples);
}

// arma::vec tn(arma::vec l, arma::vec u) {
// 
//   int n = l.n_elem;
//   arma::vec x = arma::zeros(n);
//   
//   for (int i = 0; i < n; i++) {
//     if (std::abs(u(i) - l(i)) > TOL) {
//       x(i) = trnd(l(i), u(i));
//     }
//     else if ()
//   }
//   
//   arma::uvec oidx = arma::find(arma::abs(u - l) > TOL);
//   if (oidx.n_elem > 0) {
//     x(oidx) = trnd(l(oidx), u(oidx));
//   }
//   
//   arma::uvec iidx = arma::find(arma::abs(u - l) > TOL);
// 
// }

// arma::vec trandn(arma::vec l, arma::vec u) {
//   
// 
//   
// }


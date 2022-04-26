#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "tn_moments.h"

const double NEWTON_TOL = 1e-10;
const double LOWER_BOUND_THRESH = 1e5;
const double METHOD_THRESH = 35;  // erfcinv or Newton method

// #include <vector>
// std::vector<int> primes(int n);

// arma::vec qfun(arma::vec x) {
//   arma::vec pnorm(x.n_elem);
//   for (int i = 0; i < x.n_elem; i++)
//     pnorm(i) = R::pnorm(x(i), 0., 1., false, true);
//   return arma::exp(.5 * arma::square(x) + pnorm);
// }

double qfun(double x) {
  return std::exp(.5 * std::pow(x, 2.) + R::pnorm(x, 0., 1., false, true));
}

// Newton's method for finding quantile of truncated standard normal
// arma::vec newton(arma::vec p, arma::vec l, arma::vec u) {
//   
//   arma::vec ql = qfun(l);
//   arma::vec qu(u.n_elem, arma::fill::zeros);
//   arma::uvec finite_idx = arma::find_finite(u);
//   if (finite_idx.n_elem > 0)
//     qu(finite_idx) = qfun(u(finite_idx));
//   
//   arma::vec l2 = arma::square(l);
//   arma::vec u2 = arma::square(u);
//   
//   // initial value for Newton iteration
//   arma::vec x = arma::sqrt(l - 2*arma::log(1. + p*arma::expm1(l/2. - u/2.)));
//   
//   double err = arma::datum::inf;
//   while (err > NEWTON_TOL) {
//     arma::vec del = -qfun(x) + 
//       (1-p) * arma::exp(.5*(arma::square(x) - l)) * ql +
//       p * arma::exp(.5*(arma::square(x) - u)) * qu;
//     x -= del;  // Newton step
//     err = arma::max( arma::abs(del) );
//   }
//   
//   return(x);
// }

// Newton's method for finding quantile of truncated standard normal
double newton(double p, double l, double u) {
  
  double ql = qfun(l);
  double qu = 0.;
  if (std::isfinite(u)) {
    qu = qfun(u);
  }
  
  double l2 = std::pow(l, 2.);
  double u2 = std::pow(u, 2.);
  
  // initial value for Newton iteration
  double x = std::sqrt(l - 2.*std::log(1. + p * std::expm1(.5 * (l - u))));
  
  double err = arma::datum::inf;
  while (err > NEWTON_TOL) {
    double del = -qfun(x) +
      (1-p) * std::exp(.5*(std::pow(x, 2.) - l)) * ql +
      p * std::exp(.5 * (std::pow(x, 2.) - u)) * qu;
    
    // Newton step
    x -= del;  
    err = std::abs(del);
  }
  
  return(x);
}

// normal quantile computation
// arma::vec normq(arma::vec p, arma::vec l, arma::vec u) {
//   
//   int m = l.n_elem;
//   
//   arma::vec x(l.n_elem);
//   
//   arma::uvec big_idx = arma::find(l > LOWER_BOUND_THRESH);
//   arma::uvec small_idx = arma::find(l <= LOWER_BOUND_THRESH);
//     
//   // direct computation for large x
//   if (big_idx.n_elem > 0) {
//     x(big_idx) = arma::sqrt( 
//       arma::square(l(big_idx)) - 
//       2*arma::log(
//         1. + p(big_idx) * 
//         arma::expm1(arma::square(l(big_idx)) / 2. - arma::square(u(big_idx)) / 2.)
//       )
//     );
//   }
//   
//   // Newton iterations
//   if (small_idx.n_elem > 0) {
//     x(small_idx) = newton(p(small_idx), l(small_idx), u(small_idx));
//   }
//   
//   return(x);
// }

// normal quantile computation
double normq(double p, double l, double u) {
  
  double x;
  
  if (l > LOWER_BOUND_THRESH) {
    // direct computation for large x
    x = std::sqrt(
      std::pow(l, 2.) - 2.*std::log(
          1. + p * std::expm1(.5 * (std::pow(l, 2.) + std::pow(u, 2.)))
      )
    );
  }
  else {
    x = newton(p, l, u);
  }
  
  return x;
}

// arma::vec phinv(arma::vec p, arma::vec l, arma::vec u) {
//   
//   int m = u.n_elem;
//   
//   arma::vec pl(m);
//   arma::vec x(m);
//   
//   for (int i = 0; i < m; i++) {
//     
//     // flip negative limits, use normal symmetry
//     if (u(i) < 0.) {
//       l(i) = -l(i);
//       u(i) = -u(i);
//     }
//     
//     double pl = R::pnorm(l(i), 0., 1., false, false);
//     
//     x(i) = R::qnorm(
//       pl + (R::pnorm(u(i), 0., 1., false, false) - pl) * p(i), 0., 1., false, false
//     );
//     
//     // follow through with symmetry
//     if (u(i) < 0.) {
//       x(i) = -x(i);
//     }
//   }
//   
//   return x;
// }

double phinv(double p, double l, double u) {
  double x;
  
  // flip negative limits, use normal symmetry
  if (u < 0.) {
    l = -l;
    u = -u;
  }
  
  double pl = R::pnorm(l, 0., 1., false, false);
  
  x = R::qnorm(
    pl + (R::pnorm(u, 0., 1., false, false) - pl) * p, 0., 1., false, false
  );
  
  // follow through with symmetry
  if (u < 0.) {
    x = -x;
  }
  
  return(x);
}

arma::vec norminvp(arma::vec p, arma::vec l, arma::vec u) {
  
  int m = l.n_elem;
  arma::vec x(m);
  
  // arma::vec tl, tu, tp;
  
  double tl, tu, tp;
  
  for (int i = 0; i < m; i++) {
    
    // case A: 0 < p < 1
    if (p(i) < 1. && p(i) > 0.) {
      // case 1: a < l < u
      if (l(i) > METHOD_THRESH) {
        x(i) = normq(p(i), l(i), u(i));
      }
      
      // case 2: l < u < -a
      else if (u(i) < -METHOD_THRESH) {
        x(i) = normq(1-p(i), -u(i), -l(i));
      } 
      
      // case 3: else
      else {
        x(i) = phinv(p(i), l(i), u(i));
      }  
    }
    // case B: p = 1 or p = 0
    else {
      if (std::abs(p(i) - 1.) < arma::datum::eps) {
        x(i) = u(i);
      }
      else {
        x(i) = l(i);
      }
    }
  }
  
  return x;
}

arma::ivec primes(int n) {
  std::vector<int> v;
  v.push_back(2);
  for (int i = 3; i < n; i++) {
    bool prime = true;
    for (int j=2; j*j <= i; j++)
      if (i % j == 0) {
        prime = false;
        break;
      }
      if (prime)
        v.push_back(i);
  }
  
  return v;
}

// Richtmyer sequence
//' @param dim dimension
//' @param n number of richtmyer sequence points
// [[Rcpp::export]]
arma::mat richtmyer(int dim, int n) {
  
  arma::mat x(dim, n);
  
  arma::ivec prime = primes(5 * dim * std::log( (double) dim + 1) / 4);
  arma::vec q(dim);
  for (int i = 0; i < dim; i++)
    q(i) = std::sqrt( (double) prime[i] );
  
  arma::rowvec one2N(n);
  for (int i = 0; i < n; i++)
    one2N(i) = (double) (i+1);
  arma::mat qN = q * one2N;
  
  arma::vec xr = 2. * (arma::randu(dim) - .5);
  xr.transform( [](double x) { return (x+1)*0.5; });
  
  for (int j = 0; j < n; j++) {
    x.col(j) = qN.col(j) + xr;
  }
  
  x.transform( [](double x) {return std::abs(2 * (x - (int(x)) ) - 1); });
  return x;
}


// [[Rcpp::export]]
double mvnprqmc(int n, arma::mat L, arma::vec l, arma::vec u, arma::vec mu) {
  
  int d = l.n_elem;
  
  arma::mat Z(d, n, arma::fill::zeros);
  
  arma::mat x_qmc = richtmyer(d-1, n).t();
  
  arma::vec p(n);
  arma::vec col(n, arma::fill::zeros);
  arma::vec tl;
  arma::vec tu;
  
  for (int k = 0; k < (d-1); k++) {
    
    // LZ
    if (k > 0) {
      col = ( L(k, arma::span(0, k-1)) * Z.rows(0, k-1) ).t();  
    }
    
    // truncation limits
    tl = l(k) - mu(k) - col;
    tu = u(k) - mu(k) - col;
    
    // Rcpp::Rcout << "tl size: " << tl.n_elem << std::endl;
    // Rcpp::Rcout << "x_qmc size: " << x_qmc.n_rows << ", " << x_qmc.n_cols << std::endl;
    // 
    // simulate N(mu, 1) conditional on [tl, tu] via QMC
    Z.row(k) = ( mu(k) + norminvp(x_qmc.col(k), tl, tu) ).t();
    
    // update likelihood ratio
    p += lnNpr(tl, tu) + 0.5 * std::pow(mu(k), 2.) - mu(k) * Z.row(k).t();
    
  }
  
  col = ( L.row(d-1) * Z ).t();
  tl = l(d-1) - col;
  tu = u(d-1) - col;
  p += lnNpr(tl, tu);
  
  return arma::mean(arma::exp(p));
}
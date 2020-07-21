// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
double signC(double x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  }
};

// [[Rcpp::export]]
double absC(double x) {
  if (x > 0) {
    return x;
  } else if (x == 0) {
    return 0;
  } else {
    return -x;
  }
}

// [[Rcpp::export]]
arma::vec pmax_arma(arma::vec x, double bound) {
  for (int i = 0; i < x.size(); i++) {
    if (x(i) < bound) x(i) = bound;
  }
  return x;
}

// [[Rcpp::export]]
double wcmod_basic(const arma::vec p, const arma::vec n, const double alpha) {
  arma::vec l = -log(1 - pmax_arma(p, alpha));
  arma::vec cl = cumsum(l);
  arma::vec cl2 = cl;
  int sz = cl2.size();
  cl2.resize(sz+1);
  cl2(0) = 0;
  cl2.resize(sz);
  arma::vec r = (cl2+cl)/2;
  double nt = sum(n);
  double lt = sum(l);
  double t = sqrt(nt) * ((sum(n % r)) / nt - lt / 2) / lt;
  return -log10(2 * R::pnorm(-1 * absC(t),0.0,sqrt(1.0 / 12.0),1,0)) * -1 * signC(t);
}

// [[Rcpp::export]]
arma::vec wcmod_basic2(const arma::uvec idx, const arma::mat p, const arma::mat n, const double alpha) {
  arma::vec out(p.n_rows);
  
  for (int r = 0; r < p.n_rows; ++r) {
    // Sort pval.mat and counts
    const arma::vec pval_sorted = arma::conv_to<arma::vec>::from(p.row(r))(idx);
    const arma::vec counts_sorted = arma::conv_to<arma::vec>::from(n.row(r))(idx);
    
    // Calculate
    out(r) = wcmod_basic(pval_sorted, counts_sorted, alpha);
  }
  return out;
}

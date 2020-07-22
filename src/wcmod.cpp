// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppProgress)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <progress.hpp>

using namespace RcppParallel;
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

struct wcmod : public Worker
{
  // output
  arma::mat output;
  
  // other parameters
  const arma::mat var;
  const arma::mat pval;
  const arma::mat counts;
  const double alpha;
  
  // Constructors
  wcmod(const arma::mat var, arma::mat output, const arma::mat pval, const arma::mat counts, const double alpha)
    : var(var), output(output), pval(pval), counts(counts), alpha(alpha) {}
  
  // Run function
  void operator()(size_t begin, size_t end) {
    Progress p(end*end, true);
    for (int c = begin; c < end; ++c) {
      const arma::uvec idx = arma::conv_to<arma::uvec>::from(var.col(c) - 1);
      for (int r = 0; r < pval.n_rows; ++r) {
        p.increment();
        // Sort pval.mat and counts
        const arma::vec pval_sorted = arma::conv_to<arma::vec>::from(pval.row(r))(idx);
        const arma::vec counts_sorted = arma::conv_to<arma::vec>::from(counts.row(r))(idx);
        
        // Calculate
        output(r,c) = wcmod_basic(pval_sorted, counts_sorted, alpha);
      }
    }
  }
};

// [[Rcpp::export]]
arma::mat wcmodCPP(arma::mat& var, arma::mat& pval, arma::mat& counts, double alpha) {
  // Initialize output
  arma::mat output(pval.n_rows, var.n_cols);
  
  // Create worker
  wcmod Wcmod(var, output, pval, counts, alpha);
  
  // Run
  parallelFor(0, var.n_cols, Wcmod);
  
  return Wcmod.output;
}

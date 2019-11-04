#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix baseline_cont(NumericMatrix& x, const NumericVector& y) {

  int ncolno = x.ncol();


  for(int colno = 0; colno < ncolno; colno ++) {
    x(_, colno) = x(_, colno) - y(colno);
  }
  return x;
}

// [[Rcpp::export]]

arma::cube baseline_epo(Rcpp::NumericVector& x,
                        Rcpp::NumericVector& x_m) {

  Rcpp::IntegerVector x_dims = x.attr("dim");
  Rcpp::IntegerVector xm_dims = x_m.attr("dim");
  arma::dcube y(x.begin(), x_dims[0], x_dims[1], x_dims[2], false);
  arma::mat ym(x_m.begin(), xm_dims[0], xm_dims[1], false);
  int n_chans = y.n_slices;
  int n_times = y.n_rows;

  for (int subi = 0; subi < n_chans; subi ++) {
    for (int subc = 0; subc < n_times; subc ++) {
      y.slice(subi).row(subc) -=  ym.col(subi).t();
    }
  }

  return y;
}

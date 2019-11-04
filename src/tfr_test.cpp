#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List get_list2(const Rcpp::List& x,
                     const arma::cx_mat& mors,
                     const double n_kern,
                     const std::string output) {

  Rcpp::List jj(x.length());
  Rcpp::Function rfft("mvfft");
  int npoints = mors.n_rows;
  double hmx = n_kern/2;
  int hmz = floor(hmx);
  for (int i = 0; i < x.length(); i++) {
    Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(x[i]);
    int no_rows = df.nrows();
    NumericMatrix y(no_rows, df.size());

    for (int j = 0; j < df.size() ; j++) {
      y(_,j) = NumericVector(df[j]);
    }

    arma::mat yk(y.begin(), no_rows, df.size(), false);
    yk.resize(npoints, df.size());
    ComplexMatrix yf(rfft(yk));

    arma::cx_mat yx(reinterpret_cast<std::complex<double>*>(yf.begin()), npoints, df.size());

    int timez = hmz + no_rows - 1;

    arma::cx_cube yz(no_rows, yx.n_cols,
                  mors.n_cols);

    for (arma::uword k = 0; k < yx.n_cols; k++) {
      arma::cx_mat cg(arma::ifft(mors.each_col() % yx.col(k)));
      yz.col(k) = cg.rows(hmz, timez);
    }

    if (output == "power") {
      jj[i] = arma::square(arma::abs(yz));
    } else {
      jj[i] = yz;
    }
  }

  return(jj);
}


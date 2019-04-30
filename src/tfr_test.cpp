#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List get_list2(const Rcpp::List& x,
                     const arma::cx_mat& mors,
                     const double n_kern,
                     const std::string output) {

  Rcpp::Function rfft("mvfft");
  Rcpp::List jj(x.length());
  int npoints = mors.n_rows;
  double hmx = n_kern/2;
  int hmz = floor(hmx);
  for (int i = 0; i < x.length(); i++) {
    Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(x[i]);
    int no_rows = df.nrows();
    int end_point = no_rows - 1;
    int df_cols = df.size();
    NumericMatrix y(no_rows, df_cols);

    for (int j = 0; j < df_cols ; j++) {
      y(_,j) = NumericVector(df[j]);
    }

    arma::mat yj(npoints, df_cols, arma::fill::zeros);

    for (int pj = 0; pj < df_cols ; pj++) {
      for (int iro = 0; iro < no_rows; iro++) {
        yj(iro, pj) += y(iro, pj);
      }
    }

    // Rcpp::ComplexMatrix cm = rfft(yj);
    // arma::cx_mat yx(reinterpret_cast<std::complex<double>*>(cm.begin()), cm.nrow(), cm.ncol());
    arma::cx_mat yx(arma::fft(yj));

    int time_end = yx.n_rows - hmz;
    int timez = hmz + no_rows - 1;
    int final_points = time_end - hmz;

    arma::cx_cube yz(yx.n_rows, yx.n_cols, mors.n_cols, arma::fill::zeros);
    for (arma::uword mm = 0; mm < yx.n_cols; mm++) {
      yz.col(mm) = arma::ifft(mors.each_col() % yx.col(mm));
      //yz.col(mm) = arma::ifft(yz.col(mm));
    }

    jj[i] = yz;

  }
  return(jj);
}


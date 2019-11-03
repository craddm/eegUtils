#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List get_list(const Rcpp::List& x,
                    const arma::cx_mat& mors,
                    const double n_kern,
                    const std::string output) {

  Rcpp::Function rfft("mvfft");
  Rcpp::List jj(x.length());
  unsigned int no_points = mors.n_rows;
  double hmx = n_kern/2;
  unsigned int start = floor(hmx);

  for (int i = 0; i < x.length(); i++) {
    Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(x[i]);
    unsigned int no_rows = df.nrows();
    unsigned int final_time = start + no_rows - 1;
    unsigned int no_chans = df.size();

    NumericMatrix y(no_rows, no_chans);

    for (unsigned int j = 0; j < no_chans ; j++) {
      y(_,j) = NumericVector(df[j]);
    }

    NumericMatrix yj(no_points, no_chans);

    for (unsigned int pj = 0; pj < no_rows ; pj++) {
        yj(pj,_) = y(pj,_);
    }

    Rcpp::ComplexMatrix cm = rfft(yj);
    arma::cx_mat yx(reinterpret_cast<std::complex<double>*>(cm.begin()),
                    no_points,
                    no_chans);

    arma::cx_cube yz(no_rows, no_chans,
                     mors.n_cols);

    for (arma::uword k = 0; k < no_chans ; k++) {
      Rcpp::ComplexMatrix cg = rfft(mors.each_col() % yx.col(k), true);
      arma::cx_mat yf(reinterpret_cast<std::complex<double>*>(cg.begin()), cg.nrow(), cg.ncol());
      yz.col(k) = yf.rows(start, final_time);
    }

    if (output == "power") {
      jj[i] = arma::square(arma::abs(yz / no_points));
    } else {
      jj[i] = yz / no_points;
    }
  }
  return(jj);
}

// [[Rcpp::export]]

Rcpp::List get_listal(const Rcpp::List& x,
                      const arma::cx_mat& mors,
                      const double n_kern,
                      const std::string output) {

  Rcpp::List jj(x.length());
  int no_points = mors.n_rows;
  int no_freqs = mors.n_cols;
  double hmx = n_kern/2;
  int hmz = floor(hmx);
  for (int i = 0; i < x.length(); i++) {

    Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(x[i]);

    int no_rows = df.nrows();
    int no_chans = df.size();
    int timez = hmz + no_rows - 1;

    NumericMatrix y(no_rows, no_chans);

    for (int j = 0; j < no_chans ; j++) {
      y(_,j) = NumericVector(df[j]);
    }

    arma::mat epoch_mat(y.begin(), no_rows, no_chans, false);
    epoch_mat.reshape(no_points, no_chans);

    arma::cx_mat epoch_fft(arma::fft(epoch_mat));

    arma::cx_cube yz(no_rows, no_chans,
                     no_freqs);
    arma::cx_mat cg(no_points, no_chans);

    for (arma::uword k = 0; k < no_freqs; k++) {
      cg = epoch_fft.each_col() % mors.col(k);
      cg = arma::ifft(cg);
      yz.slice(k) = cg.rows(hmz, timez);
    }

    if (output == "power") {
      jj[i] = arma::square(arma::abs(yz));
    } else {
      jj[i] = yz;
    }
  }

  return(jj);
}


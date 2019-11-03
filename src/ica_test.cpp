#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

arma::mat do_iter(arma::cube amps, int k, int N){

  int n_epochs = amps.n_slices;
  arma::mat Rxp(amps.n_rows, amps.n_rows, arma::fill::zeros);
  k += 1;
  for (int i = 0; i < n_epochs; i ++) {
    Rxp += amps.slice(i).cols(k, N - 1) * amps.slice(i).cols(0, N - k - 1).t() / (N - k);
  }
  Rxp = Rxp / n_epochs;
  Rxp = norm(Rxp, "fro") * Rxp;
  return(Rxp);
}

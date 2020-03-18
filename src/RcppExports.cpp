// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// baseline_cont
NumericMatrix baseline_cont(NumericMatrix& x, const NumericVector& y);
RcppExport SEXP _eegUtils_baseline_cont(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(baseline_cont(x, y));
    return rcpp_result_gen;
END_RCPP
}
// baseline_epo
arma::cube baseline_epo(Rcpp::NumericVector& x, Rcpp::NumericVector& x_m);
RcppExport SEXP _eegUtils_baseline_epo(SEXP xSEXP, SEXP x_mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type x_m(x_mSEXP);
    rcpp_result_gen = Rcpp::wrap(baseline_epo(x, x_m));
    return rcpp_result_gen;
END_RCPP
}
// do_iter
arma::mat do_iter(arma::cube amps, int k, int N);
RcppExport SEXP _eegUtils_do_iter(SEXP ampsSEXP, SEXP kSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type amps(ampsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(do_iter(amps, k, N));
    return rcpp_result_gen;
END_RCPP
}
// get_list
Rcpp::List get_list(const Rcpp::List& x, const arma::cx_mat& mors, const double n_kern, const std::string output);
RcppExport SEXP _eegUtils_get_list(SEXP xSEXP, SEXP morsSEXP, SEXP n_kernSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type mors(morsSEXP);
    Rcpp::traits::input_parameter< const double >::type n_kern(n_kernSEXP);
    Rcpp::traits::input_parameter< const std::string >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(get_list(x, mors, n_kern, output));
    return rcpp_result_gen;
END_RCPP
}
// get_listal
Rcpp::List get_listal(const Rcpp::List& x, const arma::cx_mat& mors, const double n_kern, const std::string output);
RcppExport SEXP _eegUtils_get_listal(SEXP xSEXP, SEXP morsSEXP, SEXP n_kernSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type mors(morsSEXP);
    Rcpp::traits::input_parameter< const double >::type n_kern(n_kernSEXP);
    Rcpp::traits::input_parameter< const std::string >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(get_listal(x, mors, n_kern, output));
    return rcpp_result_gen;
END_RCPP
}
// get_list2
Rcpp::List get_list2(const Rcpp::List& x, const arma::cx_mat& mors, const double n_kern, const std::string output);
RcppExport SEXP _eegUtils_get_list2(SEXP xSEXP, SEXP morsSEXP, SEXP n_kernSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type mors(morsSEXP);
    Rcpp::traits::input_parameter< const double >::type n_kern(n_kernSEXP);
    Rcpp::traits::input_parameter< const std::string >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(get_list2(x, mors, n_kern, output));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_eegUtils_baseline_cont", (DL_FUNC) &_eegUtils_baseline_cont, 2},
    {"_eegUtils_baseline_epo", (DL_FUNC) &_eegUtils_baseline_epo, 2},
    {"_eegUtils_do_iter", (DL_FUNC) &_eegUtils_do_iter, 3},
    {"_eegUtils_get_list", (DL_FUNC) &_eegUtils_get_list, 4},
    {"_eegUtils_get_listal", (DL_FUNC) &_eegUtils_get_listal, 4},
    {"_eegUtils_get_list2", (DL_FUNC) &_eegUtils_get_list2, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_eegUtils(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

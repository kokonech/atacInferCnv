// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// aggregate_bins_cpp_sparse
NumericMatrix aggregate_bins_cpp_sparse(IntegerVector targ_bins, IntegerVector query_hits, IntegerVector subject_hits, const Eigen::SparseMatrix<double>& signal_matrix, int num_bins);
RcppExport SEXP _atacInferCnv_aggregate_bins_cpp_sparse(SEXP targ_binsSEXP, SEXP query_hitsSEXP, SEXP subject_hitsSEXP, SEXP signal_matrixSEXP, SEXP num_binsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type targ_bins(targ_binsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type query_hits(query_hitsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type subject_hits(subject_hitsSEXP);
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type signal_matrix(signal_matrixSEXP);
    Rcpp::traits::input_parameter< int >::type num_bins(num_binsSEXP);
    rcpp_result_gen = Rcpp::wrap(aggregate_bins_cpp_sparse(targ_bins, query_hits, subject_hits, signal_matrix, num_bins));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_atacInferCnv_aggregate_bins_cpp_sparse", (DL_FUNC) &_atacInferCnv_aggregate_bins_cpp_sparse, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_atacInferCnv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

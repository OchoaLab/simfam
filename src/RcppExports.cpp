// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// draw_allele
int draw_allele(int x);
RcppExport SEXP _simfam_draw_allele(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(draw_allele(x));
    return rcpp_result_gen;
END_RCPP
}
// draw_geno_fam_cpp
IntegerMatrix draw_geno_fam_cpp(IntegerMatrix X_in, IntegerVector i_founder_in, IntegerVector i_founder_out, IntegerVector i_child, IntegerVector i_pat, IntegerVector i_mat);
RcppExport SEXP _simfam_draw_geno_fam_cpp(SEXP X_inSEXP, SEXP i_founder_inSEXP, SEXP i_founder_outSEXP, SEXP i_childSEXP, SEXP i_patSEXP, SEXP i_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type X_in(X_inSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type i_founder_in(i_founder_inSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type i_founder_out(i_founder_outSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type i_child(i_childSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type i_pat(i_patSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type i_mat(i_matSEXP);
    rcpp_result_gen = Rcpp::wrap(draw_geno_fam_cpp(X_in, i_founder_in, i_founder_out, i_child, i_pat, i_mat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_simfam_draw_allele", (DL_FUNC) &_simfam_draw_allele, 1},
    {"_simfam_draw_geno_fam_cpp", (DL_FUNC) &_simfam_draw_geno_fam_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_simfam(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

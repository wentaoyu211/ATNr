// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif


RcppExport SEXP _rcpp_module_boot_ScaledModule();
RcppExport SEXP _rcpp_module_boot_Scaled_loopsModule();
RcppExport SEXP _rcpp_module_boot_UnscaledModule();
RcppExport SEXP _rcpp_module_boot_Unscaled_loopsModule();
RcppExport SEXP _rcpp_module_boot_Unscaled_nutsModule();
RcppExport SEXP _rcpp_module_boot_Unscaled_nuts_loopsModule();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_ScaledModule", (DL_FUNC) &_rcpp_module_boot_ScaledModule, 0},
    {"_rcpp_module_boot_Scaled_loopsModule", (DL_FUNC) &_rcpp_module_boot_Scaled_loopsModule, 0},
    {"_rcpp_module_boot_UnscaledModule", (DL_FUNC) &_rcpp_module_boot_UnscaledModule, 0},
    {"_rcpp_module_boot_Unscaled_loopsModule", (DL_FUNC) &_rcpp_module_boot_Unscaled_loopsModule, 0},
    {"_rcpp_module_boot_Unscaled_nutsModule", (DL_FUNC) &_rcpp_module_boot_Unscaled_nutsModule, 0},
    {"_rcpp_module_boot_Unscaled_nuts_loopsModule", (DL_FUNC) &_rcpp_module_boot_Unscaled_nuts_loopsModule, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_ATNr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

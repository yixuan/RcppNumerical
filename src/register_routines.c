#include <cuba.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP _RcppNumerical_fastLR_(SEXP xSEXP, SEXP ySEXP, SEXP startSEXP, SEXP eps_fSEXP, SEXP eps_gSEXP, SEXP maxitSEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_RcppNumerical_fastLR_", (DL_FUNC) &_RcppNumerical_fastLR_, 6},
    {NULL, NULL, 0}
};

void R_init_RcppNumerical(DllInfo *info)
{
    R_RegisterCCallable("RcppNumerical", "Cuhre", (DL_FUNC) Cuhre);
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}

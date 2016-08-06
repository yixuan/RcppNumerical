#include <cuba.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_RcppNumerical(DllInfo *info)
{
    R_RegisterCCallable("RcppNumerical", "Cuhre", (DL_FUNC) Cuhre);
}

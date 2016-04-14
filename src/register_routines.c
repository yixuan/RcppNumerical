#include <lbfgs.h>
#include <cuba.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_RcppNumerical(DllInfo *info)
{
    R_RegisterCCallable("RcppNumerical", "lbfgs",
                        (DL_FUNC) lbfgs);
    R_RegisterCCallable("RcppNumerical", "lbfgs_parameter_init",
                        (DL_FUNC) lbfgs_parameter_init);
    R_RegisterCCallable("RcppNumerical", "Cuhre",
                        (DL_FUNC) Cuhre);
}

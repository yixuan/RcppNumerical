// Copyright (C) 2016 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPTIMIZATION_WRAPPER_H
#define OPTIMIZATION_WRAPPER_H

#include <Rcpp.h>
#include <R_ext/Rdynload.h>

#include "lbfgs.h"

namespace Numer
{

// Function lbfgs()
typedef int (*CFUN_lbfgs_TYPE)(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *ptr_fx,
    lbfgs_evaluate_t proc_evaluate,
    lbfgs_progress_t proc_progress,
    void *instance,
    lbfgs_parameter_t *param
);
// Function lbfgs_parameter_init()
typedef void (*CFUN_lbfgs_parameter_init_TYPE)(lbfgs_parameter_t *param);


inline int optim_lbfgs(
    const MFuncGrad& f, Refvec x, double& fx_opt,
    const int maxit = 300, const double& eps_f = 1e-8, const double& eps_g = 1e-6
)
{
    // Find functions
    CFUN_lbfgs_TYPE
        cfun_lbfgs = (CFUN_lbfgs_TYPE) R_GetCCallable("RcppNumerical", "lbfgs");
    CFUN_lbfgs_parameter_init_TYPE
        cfun_lbfgs_parameter_init = (CFUN_lbfgs_parameter_init_TYPE) R_GetCCallable("RcppNumerical", "lbfgs_parameter_init");

    return 0;
}


}


#endif // OPTIMIZATION_WRAPPER_H

// Copyright (C) 2016 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPTIMIZATION_WRAPPER_H
#define OPTIMIZATION_WRAPPER_H

#include <RcppEigen.h>
#include <R_ext/Rdynload.h>
#include "lbfgs.h"

namespace Numer
{


// Function type for lbfgs()
typedef int (*CFUN_lbfgs_TYPE)(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *ptr_fx,
    lbfgs_evaluate_t proc_evaluate,
    lbfgs_progress_t proc_progress,
    void *instance,
    lbfgs_parameter_t *param
);
// Function type for lbfgs_parameter_init()
typedef void (*CFUN_lbfgs_parameter_init_TYPE)(lbfgs_parameter_t *param);


// Evaluation function for lbfgs()
inline lbfgsfloatval_t lbfgs_evalfun(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{
    MFuncGrad* f = (MFuncGrad*) instance;
    const Eigen::Map<const Eigen::VectorXd> xval(x, n);
    Eigen::Map<Eigen::VectorXd> grad(g, n);

    return f->f_grad(xval, grad);
}


// [RcppNumerical API] Optimization using L-BFGS algorithm
inline int optim_lbfgs(
    MFuncGrad& f, Refvec x, double& fx_opt,
    const int maxit = 300, const double& eps_f = 1e-6, const double& eps_g = 1e-5
)
{
    // Find functions
    CFUN_lbfgs_TYPE
        cfun_lbfgs = (CFUN_lbfgs_TYPE) R_GetCCallable("RcppNumerical", "lbfgs");
    CFUN_lbfgs_parameter_init_TYPE
        cfun_lbfgs_parameter_init = (CFUN_lbfgs_parameter_init_TYPE) R_GetCCallable("RcppNumerical", "lbfgs_parameter_init");

    // Prepare parameters
    lbfgs_parameter_t param;
    cfun_lbfgs_parameter_init(&param);
    param.epsilon        = eps_g;
    param.delta          = eps_f;
    param.max_iterations = maxit;
    param.max_linesearch = 100;
    param.linesearch     = LBFGS_LINESEARCH_BACKTRACKING;

    // Call main function
    int status = cfun_lbfgs(x.size(), x.data(), &fx_opt,
                            lbfgs_evalfun, NULL, &f, &param);

    return status;
}


} // namespace Numer


#endif // OPTIMIZATION_WRAPPER_H

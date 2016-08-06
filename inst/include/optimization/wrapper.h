// Copyright (C) 2016 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPTIMIZATION_WRAPPER_H
#define OPTIMIZATION_WRAPPER_H

#include <RcppEigen.h>
#include "LBFGS.h"

namespace Numer
{


class LBFGSFun
{
private:
    MFuncGrad& f;
public:
    LBFGSFun(MFuncGrad& f_) : f(f_) {}
    inline double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
    {
        return f.f_grad(x, grad);
    }
};


// [RcppNumerical API] Optimization using L-BFGS algorithm
inline int optim_lbfgs(
    MFuncGrad& f, Refvec x, double& fx_opt,
    const int maxit = 300, const double& eps_f = 1e-6, const double& eps_g = 1e-5
)
{
    // Create functor
    LBFGSFun fun(f);

    // Prepare parameters
    LBFGSpp::LBFGSParam<double> param;
    param.epsilon        = eps_g;
    param.past           = 1;
    param.delta          = eps_f;
    param.max_iterations = maxit;
    param.max_linesearch = 100;
    param.linesearch     = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING;

    // Solver
    LBFGSpp::LBFGSSolver<double> solver(param);

    int status = 0;
    Eigen::VectorXd xx(x.size());
    xx.noalias() = x;

    try {
        solver.minimize(fun, xx, fx_opt);
    } catch(const std::exception& e) {
        status = -1;
        Rcpp::warning(e.what());
    }

    x.noalias() = xx;

    return status;
}


} // namespace Numer


#endif // OPTIMIZATION_WRAPPER_H

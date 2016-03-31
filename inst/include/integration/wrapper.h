// Copyright (C) 2016 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef INTEGRATION_WRAPPER_H
#define INTEGRATION_WRAPPER_H

#include "GaussKronrodNodesWeights.h"
#include "Integrator.h"
#include "../Func.h"

namespace Numer
{


inline double integrate(
    const Func& f, const double& lower, const double& upper,
    double& err_est, int& err_code,
    const int subdiv = 100, const double& eps_abs = 1e-8, const double& eps_rel = 1e-6,
    const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod41
)
{
    Integrator<double> intgr(subdiv);
    double res = intgr.quadratureAdaptive(f, lower, upper, eps_abs, eps_rel, rule);
    err_est = intgr.estimatedError();
    err_code = intgr.errorCode();
    return res;
}


}


#endif // INTEGRATION_WRAPPER_H

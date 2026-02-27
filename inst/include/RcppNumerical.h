// Copyright (C) 2016-2026 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef RCPPNUMERICAL_H
#define RCPPNUMERICAL_H

#include <RcppEigen.h>

// Functor definition
#include "Func.h"

// Integration
#include "integration/GaussKronrodNodesWeights.h"
#include "integration/Integrator.h"
#include "integration/cuba.h"
#include "integration/wrapper.h"

// Optimization
#include "optimization/LBFGS.h"
#include "optimization/wrapper.h"


#endif // RCPPNUMERICAL_H

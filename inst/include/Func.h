// Copyright (C) 2016 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FUNC_H
#define FUNC_H

#include <RcppEigen.h>


namespace Numer
{


// Reference to a vector
typedef Eigen::Ref<Eigen::VectorXd>             Refvec;
typedef const Eigen::Ref<const Eigen::VectorXd> Constvec;


// For 1-D numerical integration
class Func
{
public:
    virtual double operator()(const double& x) const = 0;
    virtual void eval(double* x, const int n) const
    {
        for(int i = 0; i < n; i++)
            x[i] = this->operator()(x[i]);
    }

    virtual ~Func() {}
};


// For multi-dimensional numerical integration and
// optimization that does not require gradient
class MFunc
{
public:
    virtual double operator()(Constvec& x) = 0;

    virtual ~MFunc() {}
};


// For optimization that requires gradient
class MFuncGrad
{
public:
    virtual double f_grad(Constvec& x, Refvec grad) = 0;

    virtual ~MFuncGrad() {}
};


} // namespace Numer

#endif // FUNC_H

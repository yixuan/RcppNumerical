// Copyright (C) 2016 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FUNC_H
#define FUNC_H

#include <vector>


namespace Numer
{


class Func
{
public:
    virtual double operator()(const double& x) const = 0;
    virtual void   operator()(std::vector<double>& x) const
    {
        for(std::vector<double>::size_type i = 0; i < x.size(); i++)
            x[i] = this->operator()(x[i]);
    }
};


} // namespace Numer

#endif // FUNC_H

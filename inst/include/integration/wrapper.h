// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
// Copyright (C) 2019 Ralf Stubner <ralf.stubner@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef INTEGRATION_WRAPPER_H
#define INTEGRATION_WRAPPER_H

#include "GaussKronrodNodesWeights.h"
#include "Integrator.h"
#include "cuba.h"
#include "../Func.h"

namespace Numer
{



// Internal implementation
namespace detail
{

class transform_infinite: public Func
{
private:
    const Func& func;
    const double lower;
    const double upper;
    const bool lower_finite;
    const bool upper_finite;

public:
    transform_infinite(const Func& _func, double _lower, double _upper) :
        func(_func), lower(_lower), upper(_upper),
        lower_finite(lower > -std::numeric_limits<double>::infinity()),
        upper_finite(upper < std::numeric_limits<double>::infinity())
    {}

    // Map infinite interval to (0, 1)
    double operator() (const double& t) const
    {
        const double x = (1 - t) / t;
        if (upper_finite && lower_finite)
            Rcpp::stop("At least one limit must be infinite.");
        else if (lower_finite)
            return func(lower + x) / (t * t);
        else if (upper_finite)
            return func(upper - x) / (t * t);
        else
            return (func(x) + func(-x)) / (t * t);
    }
};

} // namespace detail



//
// [RcppNumerical API] 1-D numerical integration
//
inline double integrate(
    const Func& f, const double& lower, const double& upper,
    double& err_est, int& err_code,
    const int subdiv = 100, const double& eps_abs = 1e-8, const double& eps_rel = 1e-6,
    const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod41
)
{
    // Early exit if lower and upper limits are identical
    if (upper == lower)
    {
        err_est = 0.0;
        err_code = 0;
        return 0.0;
    }

    // Finite interval
    if (std::abs(upper) < std::numeric_limits<double>::infinity() &&
        std::abs(lower) < std::numeric_limits<double>::infinity())
    {
        Integrator<double> intgr(subdiv);
        double res = intgr.quadratureAdaptive(f, lower, upper, eps_abs, eps_rel, rule);
        err_est = intgr.estimatedError();
        err_code = intgr.errorCode();
        return res;
    }

    // Infinite interval
    double sign = 1.0, lb = lower, ub = upper;
    if (ub < lb)
    {
        std::swap(ub, lb);
        sign = -1.0;
    }
    detail::transform_infinite g(f, lb, ub);

    Integrator<double> intgr(subdiv);
    double res = intgr.quadratureAdaptive(g, 0.0, 1.0, eps_abs, eps_rel, rule);
    err_est = intgr.estimatedError();
    err_code = intgr.errorCode();
    return sign * res;
}

/****************************************************************************/

// Internal implementation
namespace detail
{

// Integrate R function
class RFunc: public Func
{
private:
    Rcpp::Function fun;
    Rcpp::RObject  args;
public:
    RFunc(Rcpp::Function fun_, Rcpp::RObject args_) :
        fun(fun_),
        args(args_)
    {}

    double operator()(const double& x) const
    {
        Rcpp::NumericVector xv = Rcpp::NumericVector::create(x);
        Rcpp::NumericVector res = fun(xv, args);
        if(res.length() != 1)
            Rcpp::stop("integrand must return a vector of the same length of x");

        return Rcpp::as<double>(res);
    }

    void   operator()(double* x, const int n) const
    {
        Rcpp::NumericVector xv(n);
        std::copy(x, x + n, xv.begin());
        Rcpp::NumericVector res = fun(xv, args);
        if(res.length() != n)
            Rcpp::stop("integrand must return a vector of the same length of x");

        std::copy(res.begin(), res.end(), x);
    }
};

} // namespace detail



//
// [RcppNumerical API] 1-D numerical integration for R function
//
inline double integrate(
    Rcpp::Function f, Rcpp::RObject args, const double& lower, const double& upper,
    double& err_est, int& err_code,
    const int subdiv = 100, const double& eps_abs = 1e-8, const double& eps_rel = 1e-6,
    const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod41
)
{
    Integrator<double> intgr(subdiv);
    detail::RFunc rfun(f, args);
    double res = intgr.quadratureAdaptive(rfun, lower, upper, eps_abs, eps_rel, rule);
    err_est = intgr.estimatedError();
    err_code = intgr.errorCode();
    return res;
}

/****************************************************************************/

// Internal implementation
namespace detail
{

// Function type for Cuhre()
typedef void (*CFUN_Cuhre_TYPE)(const int ndim, const int ncomp,
    integrand_t integrand, void *userdata, const int nvec,
    const cubareal epsrel, const cubareal epsabs,
    const int flags, const int mineval, const int maxeval,
    const int key,
    const char *statefile, void *spin,
    int *nregions, int *neval, int *fail,
    cubareal integral[], cubareal err[], cubareal prob[]);

// Evaluation function for Cuhre()
inline int cuhre_integrand(const int *ndim, const cubareal x[],
                           const int *ncomp, cubareal f[], void *userdata)
{
    MFunc* func = (MFunc*) userdata;
    const Eigen::Map<const Eigen::VectorXd> xval(x, *ndim);
    *f = func->operator()(xval);

    return 0;
}

// Transform function according to integral limits
class MFuncWithBound: public MFunc
{
private:
    const double    scalefac;
    MFunc&          fun;
    Constvec&       lb;
    Eigen::VectorXd range;
    Eigen::VectorXd scalex;
public:
    MFuncWithBound(MFunc& f, Constvec& lower, Constvec& upper) :
        scalefac((upper - lower).prod()),
        fun(f), lb(lower),
        range(upper - lower), scalex(lower.size())
    {}

    inline double operator()(Constvec& x)
    {
        scalex.noalias() = lb + range.cwiseProduct(x);
        return fun(scalex);
    }

    inline double scale_factor() const { return scalefac; }

};

} // namespace detail



//
// [RcppNumerical API] Multi-dimensional integration
//
inline double integrate(
    MFunc& f, Constvec& lower, Constvec& upper,
    double& err_est, int& err_code,
    const int maxeval = 1000, const double& eps_abs = 1e-6, const double& eps_rel = 1e-6
)
{
    // Find the Cuhre() function
    detail::CFUN_Cuhre_TYPE cfun_Cuhre = (detail::CFUN_Cuhre_TYPE) R_GetCCallable("RcppNumerical", "Cuhre");

    detail::MFuncWithBound fb(f, lower, upper);
    int nregions;
    int neval;
    double integral;
    double prob;

    cfun_Cuhre(lower.size(), 1, detail::cuhre_integrand, &fb, 1,
               eps_rel, eps_abs,
               4, 1, maxeval,
               0,
               NULL, NULL,
               &nregions, &neval, &err_code, &integral, &err_est, &prob);

    integral *= fb.scale_factor();
    err_est  *= std::abs(fb.scale_factor());

    return integral;
}



}  // namespace Numer


#endif // INTEGRATION_WRAPPER_H

## Rcpp Integration for Numerical Computing Libraries

### Introduction

[Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) is a
powerful tool to write fast C++ code to speed up R programs. However,
it is not easy, or at least not straightforward to compute numerical
integration or do optimization using pure C++ code inside Rcpp.

**RcppNumerical** integrates a number of open source numerical computing
libraries into Rcpp, so that users can call convenient functions to
accomplish such tasks.

### Numerical Integration

The numerical integration code contained in **RcppNumerical** is based
on the [NumericalIntegration](https://github.com/tbs1980/NumericalIntegration)
library developed by [Sreekumar Thaithara Balan](https://github.com/tbs1980).

To compute integration of a function, first define a functor inherited from
the `Func` class:

```cpp
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
```

The first function evaluates one point at a time, and the second version
overwrites each point in the vector by the corresponding function values.
Only the second function will be used by the integration code, but usually it
is easier to implement the first one.

**RcppNumerical** provides a wrapper function for the **NumericalIntegration**
library with the following interface:

```cpp
inline double integrate(
    const Func& f, const double& lower, const double& upper,
    double& err_est, int& err_code,
    const int subdiv = 100, const double& eps_abs = 1e-8, const double& eps_rel = 1e-6,
    const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod41
)
```

- `f`: The functor of integrand.
- `lower`, `upper`: Limits of integral.
- `err_est`: Estimate of the error (output).
- `err_code`: Error code (output). See `inst/include/integration/Integrator.h`
[Line 676-704](https://github.com/yixuan/RcppNumerical/blob/master/inst/include/integration/Integrator.h#L676).
- `subdiv`: Maximum number of subintervals.
- `eps_abs`, `eps_rel`: Absolute and relative tolerance.
- `rule`: Integration rule. Possible values are
`GaussKronrod{15, 21, 31, 41, 51, 61, 71, 81, 91, 101, 121, 201}`.

See a full example below, which can be compiled using the `Rcpp::sourceCpp`
function in Rcpp.

```cpp
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>
using namespace Numer;

// Integration of Beta distribution PDF over [0.3, 0.8]

class BetaPDF: public Func
{
private:
    double a;
    double b;
public:
    BetaPDF(double a_, double b_) :
        a(a_), b(b_)
    {}

    double operator()(const double& x) const
    {
        return R::dbeta(x, a, b, 0);
    }
};


// [[Rcpp::export]]
Rcpp::List numer_test()
{
    const double a = 3, b = 10;
    const double lower = 0.3, upper = 0.8;
    const double true_val = R::pbeta(upper, a, b, 1, 0) -
                            R::pbeta(lower, a, b, 1, 0);

    BetaPDF f(a, b);
    double err_est;
    int err_code;
    const double res = integrate(f, lower, upper, err_est, err_code);
    return Rcpp::List::create(
        Rcpp::Named("true") = true_val,
        Rcpp::Named("approximate") = res,
        Rcpp::Named("error_estimate") = err_est,
        Rcpp::Named("error_code") = err_code
    );
}
```

Runing the `numer_test()` function in R gives

```r
> numer_test()
$true
[1] 0.2528108

$approximate
[1] 0.2528108

$error_estimate
[1] 2.806764e-15

$error_code
[1] 0
```

### Numerical Optimization

TODO

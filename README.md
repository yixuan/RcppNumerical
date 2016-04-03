## Rcpp Integration for Numerical Computing Libraries

### Introduction

[Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) is a
powerful tool to write fast C++ code to speed up R programs. However,
it is not easy, or at least not straightforward, to compute numerical
integration or do optimization using pure C++ code inside Rcpp.

**RcppNumerical** integrates a number of open source numerical computing
libraries into Rcpp, so that users can call convenient functions to
accomplish such tasks.

### Numerical Integration

The numerical integration code contained in **RcppNumerical** is based
on the [NumericalIntegration](https://github.com/tbs1980/NumericalIntegration)
library developed by [Sreekumar Thaithara Balan](https://github.com/tbs1980),
[Mark Sauder](https://github.com/mcsauder), and Matt Beall.

To compute integration of a function, first define a functor inherited from
the `Func` class:

```cpp
class Func
{
public:
    virtual double operator()(const double& x) const = 0;
    virtual void   operator()(double* x, const int n) const
    {
        for(int i = 0; i < n; i++)
            x[i] = this->operator()(x[i]);
    }
};
```

The first function evaluates one point at a time, and the second version
overwrites each point in the array by the corresponding function values.
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
`GaussKronrod{15, 21, 31, 41, 51, 61, 71, 81, 91, 101, 121, 201}`. Rules with
larger values have better accuracy, but may involve more function calls.

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
    BetaPDF(double a_, double b_) : a(a_), b(b_) {}

    double operator()(const double& x) const
    {
        return R::dbeta(x, a, b, 0);
    }
};

// [[Rcpp::export]]
Rcpp::List integrate_test()
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

Runing the `integrate_test()` function in R gives

```r
> integrate_test()
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

Currently **RcppNumerical** contains the L-BFGS algorithm for unconstrained
minimization problems based on the
[libLBFGS](https://github.com/chokkan/liblbfgs) library
developed by [Naoaki Okazaki](http://www.chokkan.org/).

Again, one needs to first define a functor to represent the multivariate
function to be minimized.

```cpp
class MFuncGrad
{
public:
    virtual double f_grad(Constvec& x, Refvec grad) const = 0;
};
```

Here `Constvec` represents a read-only vector and `Refvec` a writable
vector. Their definitions are

```cpp
// Reference to a vector
typedef Eigen::Ref<Eigen::VectorXd>             Refvec;
typedef const Eigen::Ref<const Eigen::VectorXd> Constvec;
```

(Basically you can treat `Refvec` as a `Eigen::VectorXd` and
`Constvec` the `const` version. Using `Eigen::Ref` is mainly to avoid
memory copy. See the explanation
[here](http://eigen.tuxfamily.org/dox/classEigen_1_1Ref.html).)

The `f_grad()` member function returns the function value on vector `x`,
and overwrites `grad` by the gradient.

The wrapper function for libLBFGS is

```cpp
inline int optim_lbfgs(
    MFuncGrad& f, Refvec x, double& fx_opt,
    const int maxit = 300, const double& eps_f = 1e-8, const double& eps_g = 1e-6
)
```

- `f`: The function to be minimized.
- `x`: In: the initial guess. Out: best value of variables found.
- `fx_opt`: Out: Function value on the output `x`.
- `maxit`: Maximum number of iterations.
- `eps_f`: Algorithm stops if `|f_{k+1} - f_k| < eps_f * |f_k|`.
- `eps_g`: Algorithm stops if `|g| < eps_g * max(1, |x|)`.

Below is an example to optimize the Rosenbrock function
`f(x1, x2) = 100 * (x2 - x1^2)^2 + (1 - x1)^2`:

```cpp
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>

using namespace Numer;

// f = 100 * (x2 - x1^2)^2 + (1 - x1)^2
// True minimum: x1 = x2 = 1
class Rosenbrock: public MFuncGrad
{
public:
    double f_grad(Constvec& x, Refvec grad) const
    {
        double t1 = x[1] - x[0] * x[0];
        double t2 = 1 - x[0];
        grad[0] = -400 * x[0] * t1 - 2 * t2;
        grad[1] = 200 * t1;
        return 100 * t1 * t1 + t2 * t2;
    }
};

// [[Rcpp::export]]
Rcpp::List optim_test()
{
    Eigen::VectorXd x(2);
    x[0] = -1.2;
    x[1] = 1;
    double fopt;
    Rosenbrock f;
    int res = optim_lbfgs(f, x, fopt);
    return Rcpp::List::create(
        Rcpp::Named("xopt") = x,
        Rcpp::Named("fopt") = fopt,
        Rcpp::Named("status") = res
    );
}
```

Calling the generated R function `optim_test()` gives

```r
> optim_test()
$xopt
[1] 1.000001 1.000001

$fopt
[1] 3.545445e-13

$status
[1] 0
```

### A More Interesting Example

It may be more meaningful to look at a real application of the **RcppNumerical**
package. Below is an example to fit logistic regression using the L-BFGS
algorithm. It also demonstrates the performance of the library.

```cpp
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>

using namespace Numer;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

class LogisticReg: public MFuncGrad
{
private:
    const MapMat X;
    const MapVec Y;
public:
    LogisticReg(const MapMat x_, const MapVec y_) : X(x_), Y(y_) {}

    double f_grad(Constvec& beta, Refvec grad) const
    {
        // Negative log likelihood
        //   sum(log(1 + exp(X * beta))) - y' * X * beta

        Eigen::VectorXd xbeta = X * beta;
        const double yxbeta = Y.dot(xbeta);
        // X * beta => exp(X * beta)
        xbeta = xbeta.array().exp();
        const double f = (xbeta.array() + 1.0).log().sum() - yxbeta;

        // Gradient
        //   X' * (p - y), p = exp(X * beta) / (1 + exp(X * beta))

        // exp(X * beta) => p
        xbeta.array() /= (xbeta.array() + 1.0);
        grad.noalias() = X.transpose() * (xbeta - Y);

        return f;
    }
};

// [[Rcpp::export]]
Rcpp::NumericVector logistic_reg(Rcpp::NumericMatrix x, Rcpp::NumericVector y)
{
    const MapMat xx = Rcpp::as<MapMat>(x);
    const MapVec yy = Rcpp::as<MapVec>(y);
    // Negative log likelihood
    LogisticReg nll(xx, yy);
    // Initial guess
    Eigen::VectorXd beta(xx.cols());
    beta.fill(0.5);

    double fopt;
    int status = optim_lbfgs(nll, beta, fopt);
    if(status != 0)
        Rcpp::stop("fail to converge");

    return Rcpp::wrap(beta);
}
```

Here is the R code to test the function:

```r
set.seed(123)
n = 5000
p = 100
x = matrix(rnorm(n * p), n)
beta = runif(p)
xb = c(x %*% beta)
p = exp(xb) / (1 + exp(xb))
y = rbinom(n, 1, p)

system.time(res1 <- glm.fit(x, y, family = binomial())$coefficient)
##  user  system elapsed
## 0.339   0.004   0.342
system.time(res2 <- logistic_reg(x, y))
##  user  system elapsed
##  0.01    0.00    0.01
max(abs(res1 - res2))
## [1] 7.862271e-09
```

It is much faster than the standard `glm.fit()` function in R! (Although
`glm.fit()` calculates some other quantities besides beta.)

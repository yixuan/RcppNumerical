#include <RcppNumerical.h>

using namespace Numer;


// [[Rcpp::export]]
Rcpp::List integrate_rcpp_(
    Rcpp::Function f, Rcpp::RObject args, double lower, double upper,
    int subdiv, double eps_abs, double eps_rel, int rule
)
{
    double err_est;
    int err_code;
    double res = integrate(f, args, lower, upper, err_est, err_code,
                           subdiv, eps_abs, eps_rel,
                           static_cast<Integrator<double>::QuadratureRule>(rule));
    return Rcpp::List::create(
        Rcpp::Named("value")     = res,
        Rcpp::Named("abs.error") = err_est,
        Rcpp::Named("status")    = err_code
    );
}

#include <RcppNumerical.h>

using namespace Numer;


// [[Rcpp::export]]
Rcpp::List integrate_rcpp_(
    Rcpp::Function f, Rcpp::RObject args, double lower, double upper,
    int subdiv, double eps_abs, double eps_rel, bool stop_on_error, int rule
)
{
    static const std::string err_msg[] = {
        "OK",
        "maximum number of subdivisions allowed has been achieved",
        "occurrence of roundoff error is detected",
        "extremely bad integrand behaviour",
        "roundoff error on extrapolation",
        "divergent integral or very slowly convergent integral",
        "input is invalid",
        "limiting number of cycles has been attained"
    };

    double err_est;
    int err_code;
    double res = integrate(f, args, lower, upper, err_est, err_code,
                           subdiv, eps_abs, eps_rel,
                           static_cast<Integrator<double>::QuadratureRule>(rule));

    if(stop_on_error)
    {
        if(err_code > 0)
            Rcpp::stop(err_msg[err_code]);
    }

    return Rcpp::List::create(
        Rcpp::Named("value")     = res,
        Rcpp::Named("abs.error") = err_est,
        Rcpp::Named("message")   = err_msg[err_code]
    );
}

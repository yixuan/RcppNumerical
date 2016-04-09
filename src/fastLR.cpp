#include <RcppNumerical.h>

using namespace Numer;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

class LogisticReg: public MFuncGrad
{
private:
    const MapMat X;
    const MapVec Y;
    const int n;
    Eigen::VectorXd xbeta;  // contains X*beta
    Eigen::VectorXd prob;   // contains log(1+exp(X*beta)) and 1/(1+exp(-X*beta))
public:
    LogisticReg(const MapMat x_, const MapVec y_) :
        X(x_),
        Y(y_),
        n(X.rows()),
        xbeta(n),
        prob(n)
    {}

    double f_grad(Constvec& beta, Refvec grad)
    {
        // Negative log likelihood
        //   sum(log(1 + exp(X * beta))) - y' * X * beta
        xbeta.noalias() = X * beta;
        const double yxbeta = Y.dot(xbeta);
        // Calculate log(1 + exp(X * beta)), avoiding overflow
        for(int i = 0; i < n; i++)
            prob[i] = R::log1pexp(xbeta[i]);
        const double f = prob.sum() - yxbeta;

        // Gradient
        //   X' * (p - y), p = exp(X * beta) / (1 + exp(X * beta))
        //                   = exp(X * beta - log(1 + exp(X * beta)))
        prob = (xbeta - prob).array().exp();
        grad.noalias() = X.transpose() * (prob - Y);

        return f;
    }

    Eigen::VectorXd current_xb() const { return xbeta; }
    Eigen::VectorXd current_p()  const { return prob; }
};

// [[Rcpp::export]]
Rcpp::List fastLR_(Rcpp::NumericMatrix x, Rcpp::NumericVector y,
                   Rcpp::NumericVector start,
                   double eps_f, double eps_g, int maxit)
{
    const MapMat xx = Rcpp::as<MapMat>(x);
    const MapVec yy = Rcpp::as<MapVec>(y);
    // Negative log likelihood
    LogisticReg nll(xx, yy);
    // Initial guess
    Rcpp::NumericVector b = Rcpp::clone(start);
    MapVec beta(b.begin(), b.length());

    double fopt;
    int status = optim_lbfgs(nll, beta, fopt, maxit, eps_f, eps_g);
    if(status < 0)
        Rcpp::warning("algorithm did not converge");

    return Rcpp::List::create(
        Rcpp::Named("coefficients")      = beta,
        Rcpp::Named("fitted.values")     = nll.current_p(),
        Rcpp::Named("linear.predictors") = nll.current_xb(),
        Rcpp::Named("loglikelihood")     = -fopt,
        Rcpp::Named("converged")         = (status >= 0)
    );
}

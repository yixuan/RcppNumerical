#include <RcppNumerical.h>

using namespace Numer;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

class LogisticReg: public MFuncGrad
{
private:
    const MapMat X;
    const MapVec Y;
    Eigen::VectorXd xbeta;  // contains X*beta and its transformation
public:
    LogisticReg(const MapMat x_, const MapVec y_) : X(x_), Y(y_) {}

    double f_grad(Constvec& beta, Refvec grad)
    {
        // Negative log likelihood
        //   sum(log(1 + exp(X * beta))) - y' * X * beta

        xbeta.noalias() = X * beta;
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

    Eigen::VectorXd current_p() const { return xbeta; }
};

// [[Rcpp::export]]
Rcpp::List logistic_reg_(Rcpp::NumericMatrix x, Rcpp::NumericVector y)
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
    if(status < 0)
        Rcpp::warning("algorithm did not converge");

    return Rcpp::List::create(
        Rcpp::Named("coefficients") = beta,
        Rcpp::Named("fitted.values") = nll.current_p(),
        Rcpp::Named("loglikelihood") = -fopt,
        Rcpp::Named("converged") = (status >= 0)
    );
}

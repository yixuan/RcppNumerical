##' Fast Logistic Regression Fitting Using L-BFGS Algorithm
##'
##' \code{fastLR()} uses the L-BFGS algorithm to efficiently fit logistic
##' regression. It is in fact an application of the C++ function
##' \code{optim_lbfgs()} provided by \pkg{RcppNumerical} to perform L-BFGS
##' optimization.
##'
##' @param x The model matrix
##' @param y The response vector
##'
##' @return \code{fastLR()} returns a list with the following components:
##' \item{coefficients}{Coefficient vector}
##' \item{fitted.values}{The fitted probability values}
##' \item{linear.predictors}{The fitted values of the linear part, i.e.,
##'                          \eqn{X\hat{\beta}}{X * beta_hat}}
##' \item{loglikelihood}{The maximized log likelihood}
##' \item{converged}{Whether the optimization algorithm has converged}
##'
##' @author Yixuan Qiu \url{http://statr.me}
##'
##' @seealso \code{\link[stats]{glm.fit}()}
##'
##' @export
##'
##' @keywords models
##' @keywords regression
##'
##' @examples
##' set.seed(123)
##' n = 1000
##' p = 100
##' x = matrix(rnorm(n * p), n)
##' beta = runif(p)
##' xb = c(x %*% beta)
##' p = 1 / (1 + exp(-xb))
##' y = rbinom(n, 1, p)
##'
##' system.time(res1 <- glm.fit(x, y, family = binomial()))
##' system.time(res2 <- fastLR(x, y))
##' max(abs(res1$coefficients - res2$coefficients))
fastLR <- function(x, y)
{
    fastLR_(x, y)
}

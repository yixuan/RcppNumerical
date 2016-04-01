integrate_rcpp <- function(f, lower, upper, ..., subdivisions = 100L,
                           rel.tol = .Machine$double.eps^0.25,
                           abs.tol = rel.tol,
                           stop.on.error = TRUE,
                           rule = 4L)
{
    f0 <- function(x, args)  do.call(f, c(list(x), args))
    args <- list(...)
    integrate_rcpp_(f0, args, lower, upper,
                    subdivisions, abs.tol, rel.tol,
                    stop.on.error, rule)
}

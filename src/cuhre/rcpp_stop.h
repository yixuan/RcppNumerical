#ifndef RCPP_STOP_H
#define RCPP_STOP_H

#ifdef __cplusplus
extern "C" {
#endif

/* Use Rcpp:stop() to handle errors */
void rcpp_stop(const char* format, ...);

#ifdef __cplusplus
}
#endif

#endif // RCPP_STOP_H

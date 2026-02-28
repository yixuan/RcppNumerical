#include <Rcpp.h>
#include "cuhre/rcpp_stop.h"

extern "C" void rcpp_stop(const char* format, ...)
{
    char buf[1024];
    va_list ap;
    va_start(ap, format);
    vsnprintf(buf, (size_t)(1024), format, ap);
    va_end(ap);
    Rcpp::stop("%s", buf);
}

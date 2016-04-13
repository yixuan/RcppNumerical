/*
	cuba.h
		Prototypes for the Cuba library
		this file is part of Cuba
		last modified 13 Mar 15 th
*/

#ifndef CUBA_H
#define CUBA_H


typedef double cubareal;

	/* integrand_t is intentionally a minimalistic integrand type.
	   It includes neither the nvec and core arguments nor the
	   extra arguments passed by Vegas/Suave (weight, iter) and
	   Divonne (phase).
	   In most cases, integrand_t is just what you want, otherwise
	   simply use an explicit typecast to integrand_t in the Cuba
	   invocation. */
typedef int (*integrand_t)(const int *ndim, const cubareal x[],
  const int *ncomp, cubareal f[], void *userdata);

typedef void (*peakfinder_t)(const int *ndim, const cubareal b[],
  int *n, cubareal x[], void *userdata);

#ifdef __cplusplus
extern "C" {
#endif


void Cuhre(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const int nvec,
  const cubareal epsrel, const cubareal epsabs,
  const int flags, const int mineval, const int maxeval,
  const int key,
  const char *statefile, void *spin,
  int *nregions, int *neval, int *fail,
  cubareal integral[], cubareal err[], cubareal prob[]);


#ifdef __cplusplus
}
#endif


#endif // CUBA_H

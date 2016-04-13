/*
	common.c
		includes most of the modules
		this file is part of Cuhre
		last modified 2 Aug 13 11 th
*/


#include "ChiSquare.c"
#include "Rule.c"

/* Use int for bool. -- Yixuan */
/* static inline bool BadDimension(cThis *t) */
static inline int BadDimension(cThis *t)
{
  /* if( t->ndim > MAXDIM ) return true; */
  if( t->ndim > MAXDIM ) return 1;
  return t->ndim < 2;
}

/* Use int for bool. -- Yixuan */
/* static inline bool BadComponent(cThis *t) */
static inline int BadComponent(cThis *t)
{
  /* if( t->ncomp > MAXCOMP ) return true; */
  if( t->ncomp > MAXCOMP ) return 1;
  return t->ncomp < 1;
}


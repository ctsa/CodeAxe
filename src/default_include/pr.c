/*
 * COPYRIGHT NOTICE, LICENSE AND DISCLAIMER
 * 
 * (c) 1994 by Peter Yianilos and Eric Sven Ristad.  All rights reserved.
 *
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation without fee for not-for-profit academic or research
 * purposes is hereby granted, provided that the above copyright notice
 * appears in all copies and that both the copyright notice and this
 * permission notice and warranty disclaimer appear in supporting
 * documentation, and that neither the authors' names nor the name of any
 * entity or institution to which they are related be used in advertising
 * or publicity pertaining to distribution of the software without
 * specific, written prior permission.
 * 
 * The authors disclaim all warranties with regard to this software,
 * including all implied warranties of merchantability and fitness.  In
 * no event shall the authors or any entities or institutions to which
 * they are related be liable for any special, indirect or consequential
 * damages or any damages whatsoever resulting from loss of use, data or
 * profits, whether in an action of contract, negligence or other
 * tortious action, arising out of or in connection with the use or
 * performance of this software.
 * 
 */

#include <assert.h>
#include <stdio.h> 
#include <stdlib.h>          /* for abort() only */
#include <math.h>
#include <limits.h>
#include <float.h>
#include "pr.h"

#define ANNOUNCE_RCSID      FALSE

#if ANNOUNCE_RCSID
static char pr_rcsid[] = "$Id: pr.c,v 1.3 1996/08/16 02:55:45 pny Exp $";
#endif

/* Note: these choices seem to have no effect on precision.  */
#define MAXIMUM_ACCURACY    FALSE
#define ADD_BEFORE_COERCE   FALSE

#if MAXIMUM_ACCURACY
 extern double expm1();
 extern double log1p();
 #define EXP(x)        (expm1(x)+1)
 #define LOG(x)        (log1p((x)-1))
#else
 #define EXP(x)        (exp(x))
 #define LOG(x)        (log(x))
#endif

#ifndef M_LN2
  #define M_LN2        6.9314718055994530942E-1
#endif
#ifndef DBL_EPSILON
  #define DBL_EPSILON  2.22044604925031310000e-16
#endif

#define NATS_POS_INFINITY  HUGE_VAL    /* nats for zero probability */

double Pr_Table[PR_PRECISION+1];  /* used by pr.bench.c so can't be static */
static double   Pr_ln2;
static BOOLEAN  Pr_Is_Initialized = FALSE;

/********************************
 *     Internal Utilities       *
 ********************************/

static void pr_initialize(void)
{
  double   x;
  unsigned u;
  int      i, j;


  Pr_Is_Initialized = TRUE;

  for (i = 0, x = 1.0; i <= PR_PRECISION; ++i, x /= PR_RADIX)
    Pr_Table[i] = x;

  Pr_ln2 = LOG(2.0);

  /* The rh term of the assert() is stored in a variable, to get around
     an apparent compiler bug in gcc version 2.6.3 (under Linux).
   */
  {
	  double ln2check;

	  ln2check = log(2.0);

	  assert(Pr_ln2 == ln2check);
  }
/*  assert(log(2.0) == LOG(2.0));*/

  i = INT_MAX;
  j = INT_MIN;
  u = (unsigned)(i - j);
  if (u != UINT_MAX)    abort();
  if (sizeof(int) < 4)  abort();
  if (PR_MANT_LOW <= 0) abort();

#ifdef FLT_RADIX
  if (FLT_RADIX != 2)   abort();
#endif

#ifdef DBL_MANT_DIG
  if (DBL_MANT_DIG != (PR_PRECISION - 1)) abort();
#endif

#if ANNOUNCE_RCSID
#ifndef NDEBUG
  fprintf(stderr, "%s\n", pr_rcsid);
#else
  fprintf(stderr, "%s\t>unsafe<\n", pr_rcsid);
#endif
#endif

}

BOOLEAN pr_is_balanced(pr_t px)
{
  return(((PR_MANT_LOW <= (px).sig) && ((px).sig < PR_MANT_HIGH)) 
	 || ((px).exp == INT_MIN));
}


void pr_fprintf(pr_t px, FILE *stream)
{
  assert(pr_is_balanced(px));
  fprintf(stream, "%fb%d", px.sig, px.exp);
}

void pr_printf(pr_t px)
{
  pr_fprintf(px, stdout);
}


/********************************
 *     Basic Type Coercion      *
 ********************************/

pr_t pr_double2pr(double x)
{
  pr_t px;

  if (!Pr_Is_Initialized) pr_initialize();

  if (x == 0) 
    {
      px.exp = INT_MIN;
      px.sig = 0;
    } 
  else 
    {
      px.sig = frexp(x, &px.exp);
    }

  assert(pr_is_balanced(px));
  return(px);
}

double pr_pr2double(pr_t px)
{
  return(ldexp(px.sig, px.exp));
}


pr_t pr_nats2pr(double x)
{
  pr_t px;

  if (!Pr_Is_Initialized) pr_initialize();

  if (x == NATS_POS_INFINITY)   /* largest codelength = smallest pr */
    {
      px = pr_double2pr(0.0);
      assert(px.exp == INT_MIN);
      assert(px.sig == 0);
    } 
  else  /* valid codelength for nonzero probability */
    {
      x = (-x) / Pr_ln2;
      px.exp = floor(x);
      px.sig = exp((x - px.exp) * Pr_ln2);
      if (px.sig >= PR_MANT_HIGH) 
	{
	  px.sig /= PR_RADIX;
	  ++px.exp;
	}
    }
  assert(pr_is_balanced(px));
  return(px);
}

double pr_pr2nats(pr_t px)
{
  double a, b;

  if (pr_is_zero(px)) return(NATS_POS_INFINITY);

  /* Note: b = log1p(px.sig - 1) does not improve accuracy */
  a = px.exp * Pr_ln2;
  b = LOG(px.sig);     
  return(-(a+b));
}

pr_t pr_bits2pr(double x)
{
  if (x == NATS_POS_INFINITY)
    return(pr_zero());
  else
    return(pr_nats2pr( x * M_LN2 ));
}


double pr_pr2bits(pr_t px)
{
  if (pr_is_zero(px))
    return(NATS_POS_INFINITY);
  else
    return(pr_pr2nats(px) / M_LN2);
}


/************************************
 * Distinguished Value Constructors *
 ************************************/

pr_t pr_zero(void)
{
  return(pr_double2pr(0.0));
}

pr_t pr_unity(void)
{
  return(pr_double2pr(1.0));
}

pr_t pr_epsilon(void)
{
  return(pr_double2pr(DBL_EPSILON));
}

/********************************
 *      Basic Predicates        *
 ********************************/

BOOLEAN pr_is_zero(pr_t px)
{
  return((px.sig == 0) && (px.exp == INT_MIN));
}

BOOLEAN pr_is_unity(pr_t px)
{
  return((px.sig == 0.5) && (px.exp == 1));
}


BOOLEAN pr_is_valid(pr_t px, pr_t epsilon)
{
  pr_t pv;
  int cmp;

  cmp = pr_exact_compare_ptr(&px,&epsilon);
  if (cmp < 0)
    {
      pv = pr_double2pr(0.0);
      pr_add_ptr(&px,&epsilon,&px);
      return(pr_exact_compare_ptr(&pv,&px) <= 0);
    }
  else if (cmp > 0)
    {
      pv = pr_double2pr(1.0);
      pr_add_ptr(&pv,&epsilon,&pv);
      return(pr_exact_compare_ptr(&px,&pv) <= 0);
    }
  else /* cmp == 0 */
    {
      return(TRUE);
    }
}


/********************************
 *      Basic Comparison        *
 ********************************/

int pr_compare(pr_t px, pr_t py, pr_t tolerance)
{
  return(pr_compare_ptr(&px,&py,&tolerance));
}

int pr_exact_compare(pr_t px, pr_t py)
{
  return(pr_exact_compare_ptr(&px,&py));
}

pr_t pr_difference(pr_t px, pr_t py)
{
  pr_t pz;
  (void) pr_difference_ptr(&px,&py,&pz);
  return(pz);
}

/********************************
 *      Basic Arithmetic        *
 ********************************/

pr_t pr_add(pr_t px, pr_t py)
{
  pr_t pz;
  pr_add_ptr(&px,&py,&pz);
  return(pz);
}

pr_t pr_multiply(pr_t px, pr_t py)
{
  pr_t pz;
  pr_multiply_ptr(&px,&py,&pz);
  return(pz);
}

pr_t pr_divide(pr_t px, pr_t py)
{
  pr_t pz;
  pr_divide_ptr(&px,&py,&pz);
  return(pz);
}

pr_t pr_power(pr_t px, double n)
{
  pr_t pz;
  pr_power_ptr(&px,n,&pz);
  return(pz);
}


/********************************
 *   Pointer-Based Comparison   *
 ********************************/

int pr_compare_ptr(pr_t *px_p, pr_t *py_p, pr_t *epsilon_p)
{
  int sign;
  pr_t magnitude;

  sign = pr_difference_ptr(px_p,py_p,&magnitude);
  if (pr_exact_compare_ptr(&magnitude,epsilon_p) <= 0) 
    return(0);
  else 
    return(sign);
}

int pr_exact_compare_ptr(pr_t *px_p, pr_t *py_p)
{
  if (px_p->exp > py_p->exp) return 1;
  if (px_p->exp < py_p->exp) return -1;
  if (px_p->sig > py_p->sig) return 1;
  if (px_p->sig < py_p->sig) return -1;
  return 0;
}

int pr_difference_ptr(pr_t *px_p, pr_t *py_p, pr_t *pz_p)
{
  int sign;
  unsigned pr_diff;

  /* store negative of smaller value in pz and add to larger value */

  if (px_p->exp > py_p->exp) 
    {
      sign = 1;
      pr_diff = (unsigned) (px_p->exp - py_p->exp);
      pz_p->exp = px_p->exp;
      if (pr_diff <= PR_PRECISION)
	pz_p->sig = px_p->sig - (py_p->sig * Pr_Table[pr_diff]);
      else
	pz_p->sig = px_p->sig;
    } 
  else if (px_p->exp < py_p->exp) 
    {
      sign = -1;
      pr_diff = (unsigned) (py_p->exp - px_p->exp);
      pz_p->exp = py_p->exp;
      if (pr_diff <= PR_PRECISION)
	pz_p->sig = py_p->sig - (px_p->sig * Pr_Table[pr_diff]);
      else
	pz_p->sig = py_p->sig;
    } 
  else if (px_p->sig > py_p->sig) 
    {
      pz_p->exp = py_p->exp;
      pz_p->sig = px_p->sig - py_p->sig;
      sign = 1;
    } 
  else if (px_p->sig < py_p->sig) 
    {
      pz_p->exp = px_p->exp;
      pz_p->sig = py_p->sig - px_p->sig;
      sign = -1;
    } 
  else 
    {
      pz_p->sig = 0;      /* pz_p->exp = INT_MIN occurs below */
      sign = 0;
    }

  /* now partially balance significand/exponent to minimize rounding error */

  if (pz_p->sig < PR_MANT_LOW) 
    {
      if (pz_p->sig == 0) 
	{
	  pz_p->exp = INT_MIN;
	} 
      else 
	{
	  while (pz_p->sig < PR_MANT_LOW)
	    {
	      assert(pz_p->exp > INT_MIN);  /* else Exponent Underflow */
	      pz_p->sig *= PR_RADIX;
	      --(pz_p->exp);
	    }
	}
    }

  assert(pr_is_balanced(*pz_p));
  return(sign);
}

/********************************
 *   Pointer-Based Arithmetic   *
 ********************************/


void pr_add_ptr(pr_t *px_p, pr_t *py_p, pr_t *pz_p)
{
  unsigned pr_diff;

  assert(Pr_Is_Initialized);

  if (px_p->exp >= py_p->exp) 
    {
      pr_diff = (unsigned) (px_p->exp - py_p->exp);
      pz_p->exp = px_p->exp;
      
      if (pr_diff <= PR_PRECISION)    /* could smaller affect larger? */
	pz_p->sig = px_p->sig + (py_p->sig * Pr_Table[pr_diff]);
      else
	pz_p->sig = px_p->sig;
    } 
  else  /* py greater than px */
    {
      pr_diff = (unsigned) (py_p->exp - px_p->exp);
      pz_p->exp = py_p->exp;
      
      if (pr_diff <= PR_PRECISION)
	pz_p->sig = py_p->sig + (px_p->sig * Pr_Table[pr_diff]);
      else
	pz_p->sig = py_p->sig;
    }

  /* now partially balance exponent and significand */
  if (pz_p->sig >= PR_MANT_HIGH) 
    {
      pz_p->sig /= PR_RADIX;
      ++(pz_p->exp);
    }

  assert(pr_is_balanced(*pz_p));
}

void pr_multiply_ptr(pr_t *px_p, pr_t *py_p, pr_t *pz_p)
{
  pz_p->sig = px_p->sig * py_p->sig;
  pz_p->exp = px_p->exp + py_p->exp;

  if (pz_p->sig < PR_MANT_LOW) 
    {
      if (pz_p->sig == 0) 
	{
	  pz_p->exp = INT_MIN;
	} 
      else 
	{
	  assert(pz_p->exp > INT_MIN);  /* else Exponent Underflow */
	  pz_p->sig *= PR_RADIX;
	  --(pz_p->exp);
	}
    }

  assert(pr_is_balanced(*pz_p));
}

void pr_divide_ptr(pr_t *px_p, pr_t *py_p, pr_t *pz_p)
{
  pz_p->sig = px_p->sig / py_p->sig;
  pz_p->exp = px_p->exp - py_p->exp;

  if (pz_p->sig >= PR_MANT_HIGH) 
    {
      pz_p->sig /= PR_RADIX;
      ++(pz_p->exp);
    } 
  else if (pz_p->sig == 0) 
    {
      pz_p->exp = INT_MIN;
    }

  assert(pr_is_balanced(*pz_p));
}


void pr_power_ptr(pr_t *px_p, double n, pr_t *pz_p)
{
  double fractional, whole, base, residue;

  fractional = modf(px_p->exp * n, &whole);
  residue = EXP(Pr_ln2 * fractional);

  base = pow(px_p->sig, n) * residue;
  pz_p->sig = frexp(base, &(pz_p->exp));

  if (pz_p->sig == 0) 
    {
      pz_p->exp = INT_MIN;
    }
  else
    {
#if ADD_BEFORE_COERCE
      whole = whole + (double)pz_p->exp;
      assert((INT_MIN < whole) && (whole < INT_MAX));
      pz_p->exp = (int)whole;
#else
      pz_p->exp = pz_p->exp + (int)whole;
#endif
    }

  assert(pr_is_balanced(*pz_p));
}



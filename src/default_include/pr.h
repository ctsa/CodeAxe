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

#ifndef PR_INCLUDED
#define PR_INCLUDED

#include <stdio.h>

#ifndef COMMON_INCLUDED
   typedef int    BOOLEAN;
#endif
#ifndef FALSE
   #define FALSE  0
#endif
#ifndef TRUE
   #define TRUE   1
#endif

#define PR_RADIX        2.0
#define PR_PRECISION    (53+1)
#define PR_MANT_LOW	(1/PR_RADIX)
#define PR_MANT_HIGH    1.0

typedef struct pr_tag
{
  double  sig;
  int     exp;
} pr_t;

/* Constructors qua Type Coercion */
pr_t	   pr_double2pr(double x);
double	   pr_pr2double(pr_t px);

pr_t	   pr_nats2pr(double x);
double	   pr_pr2nats(pr_t px);

pr_t       pr_bits2pr(double x);
double     pr_pr2bits(pr_t px);

/* Distinguished Value Constructors */
pr_t       pr_zero(void);
pr_t       pr_unity(void);
pr_t       pr_epsilon(void);

/* Predicates */
BOOLEAN    pr_is_zero(pr_t px);
BOOLEAN    pr_is_unity(pr_t px);
BOOLEAN    pr_is_valid(pr_t px, pr_t epsilon);

/* Comparison */
int	   pr_compare(pr_t px, pr_t py, pr_t tolerance);
int	   pr_exact_compare(pr_t px, pr_t py);
pr_t	   pr_difference(pr_t px, pr_t py);

/* Arithmetic Operations */
pr_t       pr_add(pr_t px, pr_t py);
pr_t 	   pr_multiply(pr_t px, pr_t py);
pr_t  	   pr_divide(pr_t px, pr_t py);
pr_t       pr_power(pr_t px, double n);

/* Pointer Comparison */
int	   pr_compare_ptr(pr_t *px_p, pr_t *py_p, pr_t *tolerance);
int	   pr_exact_compare_ptr(pr_t *px_p, pr_t *py_p);
int        pr_difference_ptr(pr_t *px_p, pr_t *py_p, pr_t *pz_p);

/* Pointer Arithmetic Operations */
void       pr_add_ptr(pr_t *px_p, pr_t *py_p, pr_t *pz_p);
void       pr_multiply_ptr(pr_t *px_p, pr_t *py_p, pr_t *pz_p);
void       pr_divide_ptr(pr_t *px_p, pr_t *py_p, pr_t *pz_p);

void       pr_power_ptr(pr_t *px_p, double n, pr_t *pz_p);

/* Input/Output */
void       pr_fprintf(pr_t px, FILE *stream);
void       pr_printf(pr_t px);

BOOLEAN    pr_is_balanced(pr_t px);

#endif	/* PR_INCLUDED */


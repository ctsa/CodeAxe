// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// CodeAxe : phylogenetic analysis and simulation tools
//
//   http://www.phrap.org
//
//

// $Id: minimize_praxis.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __MINIMIZE_PRAXIS_H
#define __MINIMIZE_PRAXIS_H


/*********************************************************************/
/* 	f u n c t i o n     p r a x i s                              */
/*                                                                   */
/* praxis is a general purpose routine for the minimization of a     */
/* function in several variables. the algorithm used is a modifi-    */
/* cation of conjugate gradient search method by powell. the changes */
/* are due to r.p. brent, who gives an algol-w program, which served */
/* as a basis for this function.                                     */
/*                                                                   */
/* references:                                                       */
/*     - powell, m.j.d., 1964. an efficient method for finding       */
/*       the minimum of a function in several variables without      */
/*       calculating derivatives, computer journal, 7, 155-162       */
/*     - brent, r.p., 1973. algorithms for minimization without      */
/*       derivatives, prentice hall, englewood cliffs.               */
/*                                                                   */
/*     problems, suggestions or improvements are always wellcome     */
/*                       karl gegenfurtner   07/08/87                */
/*                                           c - version             */
/*********************************************************************/
/*                                                                   */
/* usage: min = praxis(fun, x, n);                                   */
/*                                                                   */
/*  fun        the function to be minimized. fun is called from      */
/*             praxis with x and n as arguments                      */
/*  x          a double array containing the initial guesses for     */
/*             the minimum, which will contain the solution on       */
/*             return                                                */
/*  n          an integer specifying the number of unknown           */
/*             parameters                                            */
/*  min        praxis returns the least calculated value of fun      */
/*                                                                   */
/* some additional global variables control some more aspects of     */
/* the inner workings of praxis. setting them is optional, they      */
/* are all set to some reasonable default values given below.        */
/*                                                                   */
/*   prin      controls the printed output from the routine.         */
/*             0 -> no output                                        */
/*             1 -> print only starting and final values             */
/*             2 -> detailed map of the minimization process         */
/*             3 -> print also eigenvalues and vectors of the        */
/*                  search directions                                */
/*             the default value is 1                                */
/*  tol        is the tolerance allowed for the precision of the     */
/*             solution. praxis returns if the criterion             */
/*             2 * ||x[k]-x[k-1]|| <= sqrt(macheps) * ||x[k]|| + tol */
/*             is fulfilled more than ktm times.                     */
/*             the default value depends on the machine precision    */
/*  ktm        see just above. default is 1, and a value of 4 leads  */
/*             to a very(!) cautious stopping criterion.             */
/*  step       is a steplength parameter and should be set equal     */
/*             to the expected distance from the solution.           */
/*             exceptionally small or large values of step lead to   */
/*             slower convergence on the first few iterations        */
/*             the default value for step is 1.0                     */
/*  scbd       is a scaling parameter. 1.0 is the default and        */
/*             indicates no scaling. if the scales for the different */
/*             parameters are very different, scbd should be set to  */
/*             a value of about 10.0.                                */
/*  illc       should be set to true (1) if the problem is known to  */
/*             be ill-conditioned. the default is false (0). this    */
/*             variable is automatically set, when praxis finds      */
/*             the problem to be ill-conditioned during iterations.  */
/*  maxfun     is the maximum number of calls to fun allowed. praxis */
/*             will return after maxfun calls to fun even when the   */
/*             minimum is not yet found. the default value of 0      */
/*             indicates no limit on the number of calls.            */
/*             this return condition is only checked every n         */
/*             iterations.                                           */
/*                                                                   */
/*********************************************************************/


// C
// C     PRAXIS RETURNS THE MINIMUM OF THE FUNCTION F(X,N) OF N VARIABLES
// C     USING THE PRINCIPAL AXIS METHOD.  THE GRADIENT OF THE FUNCTION IS
// C     NOT REQUIRED.
// C
// C     FOR A DESCRIPTION OF THE ALGORITHM, SEE CHAPTER SEVEN OF
// C     "ALGORITHMS FOR FINDING ZEROS AND EXTREMA OF FUNCTIONS WITHOUT
// C     CALCULATING DERIVATIVES" BY RICHARD P BRENT.
// C
// C     THE PARAMETERS ARE:
// C     T0       IS A TOLERANCE.  PRAXIS ATTEMPTS TO RETURN PRAXIS=F(X)
// C              SUCH THAT IF X0 IS THE TRUE LOCAL MINIMUM NEAR X, THEN
// C              NORM(X-X0) < T0 + SQUAREROOT(MACHEP)*NORM(X).
// C     MACHEP   IS THE MACHINE PRECISION, THE SMALLEST NUMBER SUCH THAT
// C              1 + MACHEP > 1.  MACHEP SHOULD BE 16.**-13 (ABOUT
// C              2.22D-16) FOR REAL*8 ARITHMETIC ON THE IBM 360.
// C     H0       IS THE MAXIMUM STEP SIZE.  H0 SHOULD BE SET TO ABOUT THE
// C              MAXIMUM DISTANCE FROM THE INITIAL GUESS TO THE MINIMUM.
// C              (IF H0 IS SET TOO LARGE OR TOO SMALL, THE INITIAL RATE OF
// C              CONVERGENCE MAY BE SLOW.)
// C     N        (AT LEAST TWO) IS THE NUMBER OF VARIABLES UPON WHICH
// C              THE FUNCTION DEPENDS.
// C     PRIN     CONTROLS THE PRINTING OF INTERMEDIATE RESULTS.
// C              IF PRIN=0, NOTHING IS PRINTED.
// C              IF PRIN=1, F IS PRINTED AFTER EVERY N+1 OR N+2 LINEAR
// C              MINIMIZATIONS.  FINAL X IS PRINTED, BUT INTERMEDIATE X IS
// C              PRINTED ONLY IF N IS AT MOST 4.
// C              IF PRIN=2, THE SCALE FACTORS AND THE PRINCIPAL VALUES OF
// C              THE APPROXIMATING QUADRATIC FORM ARE ALSO PRINTED.
// C              IF PRIN=3, X IS ALSO PRINTED AFTER EVERY FEW LINEAR
// C              MINIMIZATIONS.
// C              IF PRIN=4, THE PRINCIPAL VECTORS OF THE APPROXIMATING
// C              QUADRATIC FORM ARE ALSO PRINTED.
// C     X        IS AN ARRAY CONTAINING ON ENTRY A GUESS OF THE POINT OF
// C              MINIMUM, ON RETURN THE ESTIMATED POINT OF MINIMUM.
// C     F(X,N)   IS THE FUNCTION TO BE MINIMIZED.  F SHOULD BE A REAL*8
// C              FUNCTION DECLARED EXTERNAL IN THE CALLING PROGRAM.
// C     FMIN     IS AN ESTIMATE OF THE MINIMUM, USED ONLY IN PRINTING
// C              INTERMEDIATE RESULTS.
// C     THE APPROXIMATING QUADRATIC FORM IS
// C              Q(X') = F(X,N) + (1/2) * (X'-X)-TRANSPOSE * A * (X'-X)
// C     WHERE X IS THE BEST ESTIMATE OF THE MINIMUM AND A IS
// C              INVERSE(V-TRANSPOSE) * D * INVERSE(V)
// C     (V(*,*) IS THE MATRIX OF SEARCH DIRECTIONS; D(*) IS THE ARRAY
// C     OF SECOND DIFFERENCES).  IF F HAS CONTINUOUS SECOND DERIVATIVES
// C     NEAR X0, A WILL TEND TO THE HESSIAN OF F AT X0 AS X APPROACHES X0.
// C
// C     IT IS ASSUMED THAT ON FLOATING-POINT UNDERFLOW THE RESULT IS SET
// C     TO ZERO.
// C     THE USER SHOULD OBSERVE THE COMMENT ON HEURISTIC NUMBERS AFTER
// C     THE INITIALIZATION OF MACHINE DEPENDENT NUMBERS.
// C


template <typename FloatType,
          typename MinFunc>
FloatType
minimize_praxis(FloatType* x,
                MinFunc& mf,
                const FloatType tol, // = SQREPSILON,
                const FloatType scbd = 1.0,
                const FloatType step = 1.0,
                const int ktm = 1,
                const int prin = 0,
                const int maxfun = 0,
                bool illc = false);

#include "minimize_praxis.hh"

#endif

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// CodeAxe : phylogenetic analysis and simulation tools
//
//   http://www.phrap.org
//
//

// $Id: chi_sqr.cc 957 2007-10-24 22:22:01Z ctsa $

/// \file

#include "chi_sqr.h"
#include "math_util_exception.h"

#include <cmath>

#include <algorithm>


#if 0
// f -> c++ of the  the second method bundled together with AS 245
/// WOW! the orig fortran failed spectacularly on g77, but it seems ok here, portability is suspect
///
/// replaced this with clib's lgamma
static
double
lngamma(const double z){

// c
// c       Uses Lanczos-type approximation to ln(gamma) for z > 0.
// c       Reference:
// c            Lanczos, C. 'A precision approximation of the gamma
// c                    function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
// c       Accuracy: About 14 significant digits except for small regions
// c                 in the vicinity of 1 and 2.
// c
// c       Programmer: Alan Miller
// c                   CSIRO Division of Mathematics & Statistics
// c
// c	N.B. It is assumed that the Fortran compiler supports long
// c	     variable names, including the underline character.   Some
// c	     compilers will not accept the 'implicit none' statement
// c	     below.
// c
// c       Latest revision - 17 April 1988
// c

  static const double a[] = { 0.9999999999995183, 676.5203681218835,
                              -1259.139216722289, 771.3234287757674,
                              -176.6150291498386, 12.50734324009056,
                              -0.1385710331296526, 0.9934937113930748e-05,
                              0.1659470187408462e-06};

	static const double lnsqrt2pi(0.9189385332046727);

	if(z <= 0.) throw math_util_exception("lngamma z<= 0");

	double lngamma(0.);
	double tmp(z + 7.);
	for(unsigned j(8); j>=1;--j){
	  lngamma += a[j]/tmp;
	  tmp -= 1.;
  }
	lngamma += a[0];
	lngamma = std::log(lngamma)+lnsqrt2pi-(z+6.5)+(z-0.5)*std::log(z+6.5);

  return lngamma;
}
#endif



static
double
ppnd(const double p){

// C
// C       ALGORITHM AS 111, APPL.STATIST., VOL.26, 118-121, 1977.
// C
// C       PRODUCES NORMAL DEVIATE CORRESPONDING TO LOWER TAIL AREA = P.
// C
// C	See also AS 241 which contains alternative routines accurate to
// C	about 7 and 16 decimal digits.
// C

  static const double split(0.42);

  const double q(p-0.5);

  if (std::fabs(q) <= split) {  // 0.08 <= p <= 0.92

    static const double a[] = {2.50662823884,-18.61500062529,
                               41.39119773534,-25.44106049637};
    static const double b[] = {-8.47351093090,23.08336743743,
                               -21.06224101826,3.13082909833};
    const double r(q*q);
    return q*(((a[3]*r+a[2])*r+a[1])*r+a[0])/((((b[3]*r+b[2])*r+b[1])*r+b[0])*r+1.);

  } else {

    static const double c[] = {-2.78718931138,-2.29796479134,
                               4.85014127135,2.32121276858};
    static const double d[] = {3.54388924762,1.63706781897};

    double r(std::min(p,1.-p));
    if(r<=0.) throw math_util_exception("ppnd: p range error");
    r = std::sqrt(-std::log(r));

    const double sign(q<0.?-1.:1.);
    return sign*(((c[3]*r+c[2])*r+c[1])*r+c[0])/((d[1]*r+d[0])*r+1.);
  }
}




static
double
alnorm(double x,
       bool upper){

// c
// c         Algorithm AS66 Applied Statistics (1973) vol22 no.3
// c
// c       Evaluates the tail area of the standardised normal curve
// c       from x to infinity if upper is .true. or
// c       from minus infinity to x if upper is .false.
// c

  static const double ltone(7.);
  static const double utzero(18.66);
  static const double con(1.28);

  if(x<0.){
    upper = !upper;
    x=-x;
  }

  if(! (x <= ltone || (upper && x<=utzero))) return 0.;

  const double y(0.5*x*x);

  double alnorm;

  if(x <= con){
    static const double p(0.398942280444);
    static const double q(0.39990348504);
    static const double a[] = {5.75885480458,2.62433121679,5.92885724438};
    static const double b[] = {-29.8213557807,48.6959930692};

    alnorm=0.5-x*(p-q*y/(y+a[0]+b[0]/(y+a[1]+b[1]/(y+a[2]))));

  } else {
    static const double r(0.398942280385);
    static const double c[] = {-3.8052e-8,3.98064794e-4,-0.151679116635,
                               4.8385912808,0.742380924027,3.99019417011};
    static const double d[] = {1.00000615302,1.98615381364,5.29330324926,
                               -15.1508972451,30.789933034};

    alnorm=r*std::exp(-y)/(x+c[0]+d[0]/(x+c[1]+d[1]/(x+c[2]+d[2]/(x+c[3]+d[3]/(x+c[4]+d[4]/(x+c[5]))))));
  }

  if(!upper) alnorm=1.-alnorm;
  return alnorm;
}





static
double
gammad(const double x,
       const double p){

// c
// c       algorithm as239  appl. statist. (1988) vol. 37, no. 3
// c
// c       computation of the incomplete gamma integral
// c
// c       auxiliary functions required: alogam = logarithm of the gamma
// c       function, and alnorm = algorithm as66
// c

  static const double oflo(1.e+37);
  static const double tol(1.e-14);
  static const double xbig(1.e+8);
  static const double plimit(1000.);
  static const double elimit(-88.);

  // c       check that we have valid values for x and p
  // c
  if (p <= 0. ||  x < 0.) throw math_util_exception("gammad: invalid x,p range");

  if (x == 0.) return 0.;

  // c       use a normal approximation if p > plimit
  // c
  if (p > plimit){
    const double pn1(3*std::sqrt(p)*(std::pow((x/p),-3.)+1./(9.*p)-1.));
    return alnorm(pn1,false);
  }

  // c       if x is extremely large compared to p then set gammad = 1
  // c
  if (x > xbig) { return 1.; }

  if (x <= 1. || x < p){
    // c
    // c       use pearson's series expansion.
    // c       (note that p is not large enough to force overflow in alogam).
    // c       no need to test ifault on exit since p > 0.
    // c
    double arg = p*std::log(x)-x-lgamma(p+1.);
    double gammad = 1.;

    {
      double c = 1.;
      double a = p;
      do {
        a += 1.;
        c *= x/a;
        gammad += c;
      }  while(c > tol);
    }
    arg += std::log(gammad);
    if (arg >= elimit) return std::exp(arg);
    else               return 0.;

  } else {
    // c
    // c       use a continued fraction expansion
    // c
    double arg = p*std::log(x)-x-lgamma(p);
    double a = 1. - p;
    double b = a + x + 1.;
    double c = 0.;
    double pn1 = 1.;
    double pn2 = x;
    double pn3 = x+1.;
    double pn4 = x*b;
    double gammad = pn3/pn4;
    while(true){
      a += 1.;
      b += 2.;
      c += 1.;
      const double an = a*c;
      const double pn5 = b * pn3 - an * pn1;
      const double pn6 = b * pn4 - an * pn2;
      if (std::fabs(pn6) > 0.) {
        const double rn = pn5 / pn6;
        if (std::abs(gammad - rn) <= std::min(tol,tol*rn)) break;
        gammad = rn;
      }

      pn1 = pn3;
      pn2 = pn4;
      pn3 = pn5;
      pn4 = pn6;
      if (std::fabs(pn5) >= oflo){
        // c
        // c       re-scale terms in continued fraction if terms are large
        // c
        pn1 = pn1 / oflo;
        pn2 = pn2 / oflo;
        pn3 = pn3 / oflo;
        pn4 = pn4 / oflo;
      }
    }
    arg += std::log(gammad);
    if (arg >= elimit) return 1. - std::exp(arg);
    else               return 1.;
  }
}




// AS 91 straight f to c++ translation
//
static
double
ppchi2(const double p,
       const double v,
       const double g){

// c
// c        Algorithm AS 91   Appl. Statist. (1975) Vol.24, P.35
// c
// c        To evaluate the percentage points of the chi-squared
// c        probability distribution function.
// c
// c        p must lie in the range 0.000002 to 0.999998,
// c        v must be positive,
// c        g must be supplied and should be equal to
// c          ln(gamma(v/2.0))
// c
// c     Incorporates the suggested changes in AS R85 (vol.40(1),
// c     pp.233-5, 1991) which should eliminate the need for the limited
// c     range for p above, though these limits have not been removed
// c     from the routine.
// c     If IFAULT = 4 is returned, the result is probably as accurate as
// c     the machine will allow.
// c
// c     Auxiliary routines required: PPND = AS 111 (or AS 241) and
// c     GAMMAD = AS 239.
// c

  static const unsigned maxit(20);
  static const double pmin(0.000002);
  static const double pmax(0.999998);
  static const double aa(0.6931471806);
  static const double e(0.5e-6);

  static const double c_[] = { 0., //<- added for zero indexing
                               0.01, 0.222222, 0.32, 0.4, 1.24, 2.2,
                               4.67, 6.66, 6.73, 13.32, 60.0, 70.0,
                               84.0, 105.0, 120.0, 127.0, 140.0, 175.0,
                               210.0, 252.0, 264.0, 294.0, 346.0, 420.0,
                               462.0, 606.0, 672.0, 707.0, 735.0, 889.0,
                               932.0, 966.0, 1141.0, 1182.0, 1278.0, 1740.0,
                               2520.0, 5040.0 };

  // c        test arguments and initialise
  // c
  if (p < pmin || p > pmax) throw math_util_exception("ppchi2: p out of range");
  if (v <= 0.) throw math_util_exception("ppchi2: v out of range");

  const double xx(0.5*v);
  const double c(xx-1.);

  double ch;

  // c        starting approximation for small chi-squared
  // c
  if (v < -c_[5]*std::log(p)) {
    ch = std::pow((p * xx * std::exp(g + xx * aa)),(1./xx));
    if (ch < e) return ch;

    // c        starting approximation for v less than or equal to 0.32
    // c
  } else if (v <= c_[3]) {
    ch = c_[4];
    const double a = std::log(1.-p);

    double q;
    do {
      q = ch;
      const double p1 = 1.+ch*(c_[7]+ch);
      const double p2 = ch * (c_[9] + ch * (c_[8] + ch));
      const double t = -0.5 + (c_[7] + 2. * ch) /
        p1 - (c_[9] + ch * (c_[10] + 3. * ch)) / p2;
      ch -= (1. - std::exp(a + g + 0.5 * ch + c * aa) * p2 / p1) / t;
    } while(std::fabs(q / ch - 1.) > c_[1]);

    // c        call to algorithm AS 111 - note that p has been tested above.
    // c	 AS 241 could be used as an alternative.
    // c
  } else {
    const double x(ppnd(p));

    // c        starting approximation using Wilson and Hilferty estimate
    // c
    const double p1 = c_[2] / v;
    ch = v * std::pow((x * std::sqrt(p1) + 1. - p1),3);

    // c        starting approximation for p tending to 1
    // c
    if (ch > (c_[6] * v + 6.)){
      ch = -2. * (std::log(1.-p) - c * std::log(0.5 * ch) + g);
    }
  }

  // c        call to algorithm AS 239 and calculation of seven term
  // c        Taylor series
  // c
  for(unsigned i(1);i<=maxit;++i){
    const double q = ch;
    const double p1 = 0.5 * ch;
    const double p2 = p - gammad(p1, xx);
    const double t = p2 * std::exp(xx * aa + g + p1 - c * std::log(ch));
    const double b = t / ch;
    const double a = 0.5 * t - b * c;
    const double s1 = (c_[19]+a*(c_[17]+a*(c_[14]+a*(c_[13]+a*(c_[12]+c_[11]*a)))))/c_[24];
    const double s2 = (c_[24]+a*(c_[29]+a*(c_[32]+a*(c_[33]+c_[35]*a))))/c_[37];
    const double s3 = (c_[19]+a*(c_[25]+a*(c_[28]+c_[31]*a)))/c_[37];
    const double s4 = (c_[20]+a*(c_[27]+c_[34]*a)+
                       c*(c_[22]+a*(c_[30]+c_[36]*a)))/c_[38];
    const double s5 = (c_[13]+c_[21]*a+c*(c_[18]+c_[26]*a))/c_[37];
    const double s6 = (c_[15]+c*(c_[23]+c_[16]*c))/c_[38];
    ch += t*(1.+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
    if (std::fabs(q / ch - 1.) > e) break;

    if(i==maxit) throw math_util_exception("ppchi2: max iter reached");
  }

  return ch;
}




double
chi_sqr(const double prob,
        const double df) {

  const double lng(lgamma(df/2.));
  return ppchi2(prob,df,lng);
}

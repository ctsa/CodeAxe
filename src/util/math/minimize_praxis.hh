// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// CodeAxe : phylogenetic analysis and simulation tools
//
//   http://www.phrap.org
//
//

// $Id: minimize_praxis.hh 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "array_util.h"
#include "matrix_util.h"
#include "minimize_praxis.h"
#include "minimize_praxis_minfit.h"

#include <cmath>
#include <cstdlib>

#include <limits>
#include <iomanip>
#include <iostream>


namespace PRAXIS {


// sort d and v in descending order
//
template <typename FloatType>
void
sort(FloatType* d,
     FloatType* v,
     const int n){

  for (int i(0);i<n-1;++i) {
    int k(i);
    FloatType s(d[i]);
    for (int j(i+1);j<n;++j) {
      if (d[j] > s) {
	      k = j;
	      s = d[j];
      }
    }
    if (k > i) {
      d[k] = d[i];
      d[i] = s;
      for (int j=0; j<n; j++) {
        std::swap(v[j+i*n],v[j+k*n]);
      }
    }
  }
}


template <typename FloatType>
void
matprint(const char* s,
         const FloatType * v,
         const int n,
         std::ostream& os){

  os << s << "\n";
  for (int k=0; k<n; k++){
    for (int i=0; i<n; i++){
      os << std::setw(20) << std::setprecision(10) << v[k+i*n] << " ";
    }
    os << "\n";
  }
}


template <typename FloatType>
void
vecprint(const char* s,
         const FloatType* x,
         const int n,
         std::ostream& os){

  os << s << "\n";
  for (int i=0; i<n; i++)
    os << std::setw(20) << std::setprecision(10) << x[i] << " ";
  os << "\n";
}


// print a line of traces
template <typename FloatType>
void
print(const FloatType fx,
      const int n_func_calls,
      const int n_line_searches,
      const FloatType* x,
      const int n,
      std::ostream& os){

   os << "\n";
   os << "... chi square reduced to ... "
             << std::setw(20) << std::setprecision(10) << fx << "\n";
   os << "... after " << n_func_calls << " function calls ...\n";
   os << "... including " << n_line_searches << " linear searches ...\n";
   vecprint("... current values of x ...", x, n,os);
}


// C...FLIN IS THE FUNCTION OF ONE REAL VARIABLE L THAT IS MINIMIZED
// C   BY THE SUBROUTINE MIN...
template <typename FloatType,
          typename MinFunc>
FloatType
flin(const FloatType step_size,
     const int j,
     const int n,
     const FloatType* x,
     const FloatType* v,
     const FloatType qd0,
     const FloatType qd1,
     const FloatType* q0,
     const FloatType* q1,
     int& n_func_calls,
     MinFunc& mf,
     FloatType* flin_ws){

  FloatType* tflin(flin_ws);

  if (j != -1) {		/* linear search */
    for (int i=0; i<n; i++)
      tflin[i] = x[i] + step_size *v[i+j*n];
  } else {			/* search along parabolic space curve */
    const FloatType qa = step_size*(step_size-qd1)/(qd0*(qd0+qd1));
    const FloatType qb = (step_size+qd0)*(qd1-step_size)/(qd0*qd1);
    const FloatType qc = step_size*(step_size+qd0)/(qd1*(qd0+qd1));
    for (int i=0; i<n; i++)
      tflin[i] = qa*q0[i]+qb*x[i]+qc*q1[i];
  }
  n_func_calls++;
  return mf.val(tflin);
}


// C...THE SUBROUTINE MIN MINIMIZES F FROM X IN THE DIRECTION V(*,J) UNLESS
// C   J IS LESS THAN 1, WHEN A QUADRATIC SEARCH IS MADE IN THE PLANE
// C   DEFINED BY Q0,Q1,X.
// C   D2 IS EITHER ZERO OR AN APPROXIMATION TO HALF F".
// C   ON ENTRY, X1 IS AN ESTIMATE OF THE DISTANCE FROM X TO THE MINIMUM
// C   ALONG V(*,J) (OR, IF J=0, A CURVE).  ON RETURN, X1 IS THE DISTANCE
// C   FOUND.
// C   IF FK=.TRUE., THEN F1 IS FLIN(X1).  OTHERWISE X1 AND F1 ARE IGNORED
// C   ON ENTRY UNLESS FINAL FX IS GREATER THAN F1.
// C   NITS CONTROLS THE NUMBER OF TIMES AN ATTEMPT WILL BE MADE TO HALVE
// C   THE INTERVAL
template <typename FloatType,
          typename MinFunc>
void
min(const int j, const int nits, FloatType& d2, FloatType& step_size, FloatType f1, const bool fk,
    const int n, FloatType& fx,
    const FloatType macheps,const FloatType m4, const FloatType dmin,
    const FloatType m2, const FloatType ldt, const FloatType small,
    const FloatType h, int& n_func_calls, MinFunc& mf, int& n_line_searches,
    const FloatType tolx,
    FloatType* x, const FloatType* v,
    const FloatType qd0, const FloatType qd1,
    const FloatType* q0, const FloatType* q1,
    FloatType* tmp_ws){

  const FloatType start_f1 = f1;
  const FloatType start_step_size = step_size;
  int k = 0;
  FloatType xm = 0.0;
  FloatType f_min = fx;
  FloatType f0 = fx;
  bool dz = d2 < macheps;


  /* find step size */
  FloatType t2;
  {
    FloatType s(array_norm2(x,n));
    if (dz) {
      t2 = m4*std::sqrt(std::fabs(fx)/dmin + s*ldt) + m2*ldt;
    } else {
      t2 = m4*std::sqrt(std::fabs(fx)/(d2) + s*ldt) + m2*ldt;
    }
    s = s*m4 + tolx;
    if (dz && t2 > s) t2 = s;
  }
  if (t2 < small) t2 = small;
  if (t2 > 0.01*h) t2 = 0.01 * h;

  if (fk && f1 <= f_min) {
    xm = step_size;
    f_min = f1;
  }
  if (!fk || std::fabs(step_size) < t2) {
    step_size = (step_size > 0 ? t2 : -t2);
    f1 = flin(step_size, j,n,x,v,qd0,qd1,q0,q1,n_func_calls,mf,tmp_ws);
  }
  if (f1 <= f_min) {
    xm = step_size;
    f_min = f1;
  }

  FloatType x2,f2;
nnext:
  if (dz) {
    //C...EVALUATE FLIN AT ANOTHER POINT AND ESTIMATE THE SECOND DERIVATIVE...
    x2 = (f0 < f1 ? -(step_size) : 2*(step_size));
    f2 = flin(x2, j, n,x,v,qd0,qd1,q0,q1,n_func_calls,mf,tmp_ws);
    if (f2 <= f_min) {
      xm = x2;
      f_min = f2;
    }
    d2 = (x2*(f1-f0) - (step_size)*(f2-f0))/((step_size)*x2*((step_size)-x2));
  }
  {
    //C...ESTIMATE THE FIRST DERIVATIVE AT 0...
    const FloatType d1 = (f1-f0)/(step_size) - step_size*d2;
    dz = true;
    //C...PREDICT THE MINIMUM...
    if (d2 <= small) {
      x2 = (d1 < 0 ? h : -h);
    } else {
      x2 = - 0.5*d1/(d2);
    }
  }
  if (std::fabs(x2) > h) {
    x2 = (x2 > 0 ? h : -h);
  }

  while(true){
    //C...EVALUATE F AT THE PREDICTED MINIMUM...
    f2 = flin(x2, j, n, x,v,qd0,qd1,q0,q1, n_func_calls,mf,tmp_ws);
    if ((k >= nits) || (f2 <= f0)) break;

    //C...NO SUCCESS, SO TRY AGAIN...
    k++;
    if ((f0 < f1) && (step_size*x2 > 0.0)) goto nnext;
    x2 *= 0.5;
  }
  n_line_searches++;
  if (f2 > f_min) x2 = xm;
  else         f_min = f2;

  //C...GET A NEW ESTIMATE OF THE SECOND DERIVATIVE...
  if (fabs(x2*(x2-step_size)) > small) {
    d2 = (x2*(f1-f0) - step_size*(f_min-f0))/(step_size*x2*(step_size-x2));
  } else {
    if (k > 0) d2 = 0;
  }
  if (d2 <= small) d2 = small;
  step_size = x2;
  fx = f_min;
  if (start_f1 < fx) {
    fx = start_f1;
    step_size = start_step_size;
  }

  //C...UPDATE X FOR LINEAR BUT NOT PARABOLIC SEARCH...
  if (j != -1) {
    for (int i=0; i<n; i++) {
      x[i] += (step_size)*v[i+j*n];
    }
  }
}


// look for a minmum along the curve q0, q1, x
template <typename FloatType,
          typename MinFunc>
void
quad_min(FloatType& qf1, FloatType& qd1, FloatType* x, FloatType* q1,
         const int n, FloatType& fx,
         const FloatType macheps,const FloatType m4, const FloatType dmin,
         const FloatType m2, const FloatType ldt, const FloatType small,
         const FloatType h, int& n_func_calls, MinFunc& mf, int& n_line_searches,
         const FloatType tolx,const FloatType* v,FloatType& qd0,FloatType* q0,
         FloatType* tmp_ws){

  FloatType l;
  FloatType s = fx;
  fx = qf1;
  qf1 = s;
  qd1 = 0.0;
  for (int i=0; i<n; i++) {
    s = x[i];
    l = q1[i];
    x[i] = l;
    q1[i] = s;
    qd1 = qd1 + (s-l)*(s-l);
  }
  s = 0.0;
  qd1 = std::sqrt(qd1);
  l = qd1;

  FloatType qa,qb,qc;
  if (qd0>0.0 && qd1>0.0 && n_line_searches>=3*n*n) {
    min(-1, 2, s, l, qf1,true,n,fx,macheps,m4,dmin,m2,ldt,small,h,n_func_calls,mf,n_line_searches, tolx,x,v,qd0,qd1,q0,q1,tmp_ws);
    qa = l*(l-qd1)/(qd0*(qd0+qd1));
    qb = (l+qd0)*(qd1-l)/(qd0*qd1);
    qc = l*(l+qd0)/(qd1*(qd0+qd1));
  }
  else {
    fx = qf1;
    qa = qb = 0.0;
    qc = 1.0;
  }
  qd0 = qd1;
  for (int i=0; i<n; i++) {
    s = q0[i];
    q0[i] = x[i];
    x[i] = qa*s + qb*x[i] + qc*q1[i];
  }
}

} // namespace PRAXIS

template <typename FloatType,
          typename MinFunc>
FloatType
minimize_praxis(FloatType* x,
                MinFunc& mf,
                const FloatType tol,
                const FloatType scbd,
                const FloatType step,
                const int ktm,
                const int prin,
                const int maxfun,
                bool is_ill_conditioned){

  using namespace PRAXIS;

  static const int min_step_factor(1);

  static const FloatType macheps = std::numeric_limits<FloatType>::epsilon();
  static const FloatType small = macheps*macheps;
  static const FloatType vsmall = small*small;
  static const FloatType large = 1.0/small;
  static const FloatType vlarge = 1.0/vsmall;
  static const FloatType m2 = std::sqrt(macheps);
  static const FloatType m4 = std::sqrt(m2);

  static std::ostream& log_os(std::cerr);

  const int n(mf.dim());

  const unsigned ws_size(sizeof(FloatType)*n*(7+n));
  char* workspace(new char[ws_size]);
  FloatType* d(reinterpret_cast<FloatType*>(workspace));
  FloatType* vec_start(d+n);
  FloatType* vec_newdir(vec_start+n);
  FloatType* z(vec_newdir+n);
  FloatType* q0(z+n);
  FloatType* q1(q0+n);
  FloatType* tmp_ws(q1+n);
  FloatType* v(tmp_ws+n);

  FloatType ldfac = (is_ill_conditioned ? 0.1 : 0.01);
  int n_line_searches = 0;
  int kt = 0;
  int n_func_calls = 1;
  FloatType fx = mf.val(x);
  FloatType qf1 = fx;
  FloatType dmin = small;
  FloatType t2 = small + std::fabs(tol);
  const FloatType tolx = t2;
  FloatType h = step;


  if (h < 100.0*tolx) h = 100.0*tolx;
  FloatType ldt = h;

  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
       v[i+j*n] = (i == j ? 1.0 : 0.0);
    }
  }

  d[0] = 0.0;
  FloatType qd0 = 0.0;
  for (int i=0; i<n; i++) q1[i] = x[i];

  if (prin > 1) {
    log_os << "\n------------- enter function praxis -----------\n";
    log_os << "... current parameter settings ...\n";
    log_os << "... scaling ... " << std::setw(20) << std::setprecision(10) << scbd << "\n";
    log_os << "...   tol   ... " << std::setw(20) << std::setprecision(10) << tolx << "\n";
    log_os << "... maxstep ... " << std::setw(20) << std::setprecision(10) << h << "\n";
    log_os << "...   is_ill_conditioned  ... " << std::setw(20) << is_ill_conditioned << "\n";
    log_os << "...   ktm   ... " << std::setw(20) << ktm << "\n";
    log_os << "... maxfun  ... " << std::setw(20) << maxfun << "\n";
  }
  if (prin) print(fx,n_func_calls,n_line_searches,x,n,log_os);

  int kl;
  FloatType sf,s,qd1,lds,dn;

  //C .....THE MAIN LOOP STARTS HERE.....
  while(true){
    sf = d[0];
    s = d[0] = 0.0;

    /* minimize along first direction */
    min(0, 2*min_step_factor, d[0], s, fx, false, n, fx, macheps,m4,dmin,m2,ldt,small,h,n_func_calls,mf,n_line_searches, tolx,x,v,qd0,qd1,q0,q1,tmp_ws);
    if (s <= 0.0)
      for (int i=0; i < n; i++)
        v[i] = -v[i];
    if ((sf <= (0.9 * d[0])) || ((0.9 * sf) >= d[0]))
      for (int i=1; i<n; i++)
        d[i] = 0.0;

    //C.....THE INNER LOOP STARTS HERE.....
    FloatType df,f1;
    for (int k=1; k<n; k++) {
      for (int i=0; i<n; i++) vec_start[i] = x[i];
      sf = fx;
      if(kt > 0) is_ill_conditioned=true;

 nnext:
      kl = k;
      df = 0.0;
      if (is_ill_conditioned) {        /* random step to get off resolution valley */
        for (int i=0; i<n; i++) {
          z[i] = (0.1 * ldt + t2 * std::pow(10.0,static_cast<FloatType>(kt))) * (drand48() - 0.5);
          s = z[i];
          for (int j=0; j < n; j++)
            x[j] += s * v[j+i*n];
        }
        fx = mf.val(x);
        n_func_calls++;
      }
      /* minimize along non-conjugate directions */
      for (int k2=k; k2<n; k2++) {
        const FloatType sl = fx;
        s = 0.0;
        min(k2, 2*min_step_factor, d[k2], s, fx, false, n, fx, macheps, m4,dmin,m2,ldt,small,h,n_func_calls,mf,n_line_searches, tolx,x,v,qd0,qd1,q0,q1,tmp_ws);
        if (is_ill_conditioned) {
          FloatType szk = s + z[k2];
          s = d[k2] * szk*szk;
        } else {
          s = sl - fx;
        }
        if (df < s) {
          df = s;
          kl = k2;
        }
      }
      if (!is_ill_conditioned && (df < std::fabs(100. * macheps * fx))) {

        //C.....IF THERE WAS NOT MUCH IMPROVEMENT ON THE FIRST TRY, SET
        //C     IS_ILL_CONDITIONED=TRUE AND START THE INNER LOOP AGAIN.....
        //C
        is_ill_conditioned = true;
        goto nnext;
      }
      if ((k == 1) && (prin > 1))
        vecprint("\n... New Direction ...",d,n,log_os);

      /* minimize along conjugate directions */
      for (int k2=0; k2<=k-1; k2++) {
        s = 0.0;
        min(k2, 2*min_step_factor, d[k2], s, fx, false, n, fx, macheps, m4, dmin, m2,ldt,small,h,n_func_calls,mf,n_line_searches, tolx,x,v,qd0,qd1,q0,q1,tmp_ws);
      }
      f1 = fx;
      fx = sf;
      for (int i=0; i<n; i++) {
        const FloatType sl = x[i]-vec_start[i];
        x[i] = vec_start[i];
        vec_newdir[i] = sl;
      }
      lds = array_norm2(vec_newdir,n);
      if (lds > small) {

        // C
        // C.....DISCARD DIRECTION V(*,KL).
        // C     IF NO RANDOM STEP WAS TAKEN, V(*,KL) IS THE "NON-CONJUGATE"
        // C     DIRECTION ALONG WHICH THE GREATEST IMPROVEMENT WAS MADE.....
        // C
        for (int i=kl-1; i>=k; i--) {
          for (int j=0; j < n; j++)
            v[j+(i+1)*n] = v[j+i*n];
          d[i+1] = d[i];
        }
        d[k] = 0.0;
        for (int i=0; i < n; i++)
          v[i+k*n] = vec_newdir[i] / lds;

        // C
        // C.....MINIMIZE ALONG THE NEW "CONJUGATE" DIRECTION V(*,K), WHICH IS
        // C     THE NORMALIZED VECTOR:  (NEW X) - (0LD X).....
        // C
        min(k, 4*min_step_factor, d[k], lds, f1, true, n, fx, macheps, m4, dmin, m2, ldt,small,h,n_func_calls,mf,n_line_searches, tolx,x,v,qd0,qd1,q0,q1,tmp_ws);
        if (lds <= 0.0) {
          lds = -lds;
          for (int i=0; i<n; i++)
            v[i+k*n] = -v[i+k*n];
        }
      }
      ldt *= ldfac;
      if (ldt < lds)
        ldt = lds;
      if (prin > 1)
        print(fx,n_func_calls,n_line_searches,x,n,log_os);
      t2 = array_norm2(x,n);
      t2 = m2 * t2 + tolx;

      // C
      // C.....SEE WHETHER THE LENGTH OF THE STEP TAKEN SINCE STARTING THE
      // C     INNER LOOP EXCEEDS HALF THE TOLERANCE.....
      // C
      if (ldt > (0.5 * t2))
        kt = 0;
      else
        kt++;
      if (kt > ktm) goto fret;
    }

    // C
    // C     TRY QUADRATIC EXTRAPOLATION IN CASE WE ARE IN A CURVED VALLEY.
    // C
    quad_min(qf1, qd1, x, q1, n, fx, macheps, m4, dmin, m2, ldt,small,h,n_func_calls,mf,n_line_searches, tolx,v,qd0,q0,tmp_ws);
    dn = 0.0;
    for (int i=0; i<n; i++) {
      d[i] = 1./std::sqrt(d[i]);
      if (dn < d[i]) dn = d[i];
    }
    if (prin > 2)
      matprint("\n... New Matrix of Directions ...",v,n,log_os);
    for (int j=0; j<n; j++) {
      s = d[j] / dn;
      for (int i=0; i < n; i++) v[i+j*n] *= s;
    }

    // C
    // C.....SCALE THE AXES TO TRY TO REDUCE THE CONDITION NUMBER.....
    // C
    if (scbd > 1.0) {
      s = vlarge;
      for (int i=0; i<n; i++) {
        FloatType sl = 0.0;
        for (int j=0; j < n; j++)
          sl += v[i+j*n]*v[i+j*n];
        z[i] = std::sqrt(sl);
        if (z[i] < m4) z[i] = m4;
        if (s > z[i]) s = z[i];
      }
      for (int i=0; i<n; i++) {
        FloatType sl = s / z[i];
        z[i] = 1.0 / sl;
        if (z[i] > scbd) {
          sl = 1.0 / scbd;
          z[i] = scbd;
        }
        for(int j=0; j<n; ++j) v[i+j*n] *= sl;
      }
    }

    // C
    // C.....CALCULATE A NEW SET OF ORTHOGONAL DIRECTIONS BEFORE REPEATING
    // C     THE MAIN LOOP.
    // C     FIRST TRANSPOSE V FOR MINFIT:
    // C
    matrix_transpose_inplace(v,n);

    // C
    // C.....CALL MINFIT TO FIND THE SINGULAR VALUE DECOMPOSITION OF V.
    // C     THIS GIVES THE PRINCIPAL VALUES AND PRINCIPAL DIRECTIONS OF THE
    // C     APPROXIMATING QUADRATIC FORM WITHOUT SQUARING THE CONDITION
    // C     NUMBER.....
    // C
    minfit(n, macheps, vsmall, v, d, tmp_ws);

    // C
    // C.....UNSCALE THE AXES.....
    // C
    if (scbd > 1.0) {
      for (int i=0; i<n; i++) {
        s = z[i];
        for (int j=0; j<n; j++) v[i+j*n] *= s;
      }
      for (int i=0; i<n; i++) {
        s = 0.0;
        for (int j=0; j<n; j++) s += v[j+i*n]*v[j+i*n];
        s = std::sqrt(s);
        d[i] *= s;
        s = 1.0 / s;
        for (int j=0; j<n; j++) v[j+i*n] *= s;
      }
    }

    for (int i=0; i<n; i++) {
      if ((dn * d[i]) > large)
        d[i] = vsmall;
      else if ((dn * d[i]) < small)
        d[i] = vlarge;
      else
        d[i] = std::pow(dn * d[i],-2.0);
    }

    // C
    // C.....SORT THE EIGENVALUES AND EIGENVECTORS.....
    // C
    sort(d,v,n);
    dmin = d[n-1];
    if (dmin < small)
      dmin = small;
    is_ill_conditioned = (m2 * d[0]) > dmin;
    if ((prin > 2) && (scbd > 1.0))
      vecprint("\n... Scale Factors ...",z,n,log_os);
    if (prin > 2)
      vecprint("\n... Eigenvalues of A ...",d,n,log_os);
    if (prin > 2)
      matprint("\n... Eigenvectors of A ...",v,n,log_os);

    if ((maxfun > 0) && (n_func_calls > maxfun)) {
      if (prin)
        log_os << "\n... maximum number of function calls reached ...\n";
      break;
    }
  } 	 /* main loop */

fret:
  if (prin > 0) {
    vecprint("\n... Final solution is ...", x, n,log_os);
    log_os << "\n... Func value reduced to " << std::setw(20) << std::setprecision(10) << fx << " ...\n";
    log_os << "... after " << std::setw(20) << n_func_calls << " function calls.\n";
  }

  delete [] workspace;

  return(fx);
}

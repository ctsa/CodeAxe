// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// CodeAxe : phylogenetic analysis and simulation tools
//
//   http://www.phrap.org
//
//
// Copyright 2007 Christopher T Saunders (ctsa@u.washington.edu)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
//
//

// $Id: matrix_exp_diag.cc 1148 2008-03-11 00:49:49Z ctsa $

/// \file

#include "matrix_exp_diag.h"
#include "matrix_exp_diag_float_thresh.h"
#include "matrix_invert.h"
#include "matrix_util.h"
#include "math_util_exception.h"
#include "../general/die.h"
#include "../general/log.h"

#include <cstring>

#include <ostream>
#include <limits>


#ifdef USE_LAPACK
#include "clapack_util.h"


// get eigenvectors/values of in_mat
//
// ** in_mat is destroyed **
//
void
matrix_exp_diag_get_evec(double* eval_real,  // dim
                         double* eval_cplx,  // dim
                         double* evec,       // dim*dim
                         double* in_mat,     // dim*dim
                         const unsigned N){

  using namespace CLAPACK;

  integer n = static_cast<integer>(N);

  double* dummy=0;
  integer idummy=1;

  double dummy_work = 0;
  integer work_size = -1;

  integer info = 0;

  // test in_mat for symmetry
  if( matrix_is_symm(in_mat,N) ){

    dsyev_("V","U", &n, in_mat, &n, eval_real,&dummy_work,&work_size,&info);

    if( info != 0 ){
      die("dsyev (symmetric eigensolver) worksize failure in matrix_exp_diag_get_evec");
    }

    work_size = static_cast<integer>(dummy_work);
    double* work = new double[work_size];

    dsyev_("V","U", &n, in_mat, &n, eval_real,work,&work_size,&info);

    if( info != 0 ){
      die("dsyev (symmetric eigensolver) failure in matrix_exp_diag_get_evec");
    }

    delete[] work;

    array_zero(eval_cplx,N);
    for(unsigned i(0);i<N*N;++i) evec[i] = in_mat[i];

  } else {
    // query workspace size
#ifndef BALANCED_GEEV
    dgeev_("N","V", &n, in_mat, &n, eval_real, eval_cplx,
           dummy, &idummy, evec, &n,&dummy_work,&work_size,&info);
#else
    integer ilo,ihi;
    double* scale(new double[n]);
    double* rconde(new double[n]);
    double* rcondv(new double[n]);
    double abnrm;
    integer* iwork;
    dgeevx_("B","N","V","N", &n, in_mat, &n, eval_real, eval_cplx,
            dummy, &idummy, evec, &n, &ilo, &ihi, scale, &abnrm,
            rconde, rcondv, &dummy_work,&work_size,iwork,&info);
#endif

    if( info != 0 ){
      die("dgeev (general eigensolver) worksize failure in matrix_exp_diag_get_evec");
    }

    work_size = static_cast<integer>(dummy_work);
    double* work = new double[work_size];

#ifndef BALANCED_GEEV
    // get eigenvectors
    dgeev_("N","V", &n, in_mat, &n, eval_real, eval_cplx,
           dummy, &idummy, evec, &n,work,&work_size,&info);
#else
    dgeevx_("B","N","V","N", &n, in_mat, &n, eval_real, eval_cplx,
            dummy, &idummy, evec, &n, &ilo, &ihi, scale, &abnrm,
            rconde, rcondv, work,&work_size,iwork,&info);
#endif

    if( info != 0 ){
      die("dgeev (general eigensolver) failure in matrix_exp_diag_get_evec");
    }

    delete[] work;
  }
}



void
matrix_exp_diag_prep(matrix_exp_diag_prepdata& pd,
                     const float* in_mat){

  const unsigned N = pd.dim;
  double* in_mat_d = new double[N*N];

  for(unsigned i(0);i<N*N;++i) { in_mat_d[i] = in_mat[i]; }

  matrix_exp_diag_prep(pd,in_mat_d);

  delete [] in_mat_d;
}

#endif

void
matrix_exp_diag_scale(float* out_mat,     // dim*dim
                      const float scale,
                      matrix_exp_diag_prepdata& pd,
                      const bool is_get_deriv){

  const unsigned N = pd.dim;
  double* out_mat_d = new double[N*N];

  for(unsigned i(0);i<N*N;++i) out_mat_d[i] = out_mat[i];
  matrix_exp_diag_scale(out_mat_d,static_cast<double>(scale),pd,is_get_deriv);
  for(unsigned i(0);i<N*N;++i) out_mat[i] = out_mat_d[i];

  delete [] out_mat_d;
}

#ifdef USE_LAPACK
void
matrix_exp_diag_prep(matrix_exp_diag_prepdata& pd,
                     const double* in_mat){  // dim*dim

  static const double dep(std::numeric_limits<double>::epsilon());

  using namespace CLAPACK;

  const unsigned N = pd.dim;
  const unsigned N2 = N*N;

  // get rate matrix eigenvectors::
  //

  // make copy because lapack destroys input matrix
  double* in_mat_copy = new double[N2];
  memcpy(in_mat_copy,in_mat,N2*sizeof(double));

  double* eval_real = pd.eval;
  double* eval_cplx = new double[N];
  matrix_exp_diag_get_evec(eval_real,eval_cplx,pd.evec,in_mat_copy,N);

  delete [] in_mat_copy;

//   for(unsigned i(0);i<N;++i){
//     log_os << "eigenvals r/c: "
//               << pd.eval_real[i] << " " << pd.eval_cplx[i] << "\n";
//   }

//   log_os << "eigenvectors:\n";
//   for(unsigned i(0);i<N;++i){
//     log_os << "col "<< i << "\n";
//     for(unsigned j(0);j<N;++j){ log_os << pd.evec[j+i*N] << " ";}
//     log_os << "\n";
//   }


  pd.is_complex = false;

  // look for complex eigenvalues -- if found switch to complex matrices
  for(unsigned i(0);i<N;++i){
    const double aer(std::fabs(eval_real[i]));
    const double aec(std::fabs(eval_cplx[i]));
    if( (aer > 0. && aec > aer*dep) ||
        (aer <= 0. && aec > aer)){
      pd.is_complex = true;
      break;
    }
  }

  if( ! pd.is_complex ){
    // invert rate matrix eigenvectors::
    memcpy(pd.inv_evec,pd.evec,N2*sizeof(double));
    matrix_invert<double>(pd.inv_evec,N);
  } else {
    // decode evec into complex form:
    //
    unsigned imag_count = 0;
    for(unsigned i(0);i<N;++i){
      pd.Ceval[i] = std::complex<double>(eval_real[i],eval_cplx[i]);
      const double aer(std::fabs(eval_real[i]));
      const double aec(std::fabs(eval_cplx[i]));
      if( (aer > 0. && aec > aer*dep) ||
          (aer <= 0. && aec > aer)){
        const unsigned iodd(imag_count%2);
        if(iodd==0){
          const unsigned realindex((i-(iodd))*N);
          const unsigned imagindex(realindex+N);
          for(unsigned j(0);j<N;++j){
            pd.Cevec[j+i*N] =
              std::complex<double>(pd.evec[j+realindex],pd.evec[j+imagindex]);
          }
        } else {
          for(unsigned j(0);j<N;++j){
            pd.Cevec[j+i*N] = std::conj(pd.Cevec[j+(i-1)*N]);
          }
        }
        imag_count++;
      } else {
        for(unsigned j(0);j<N;++j){
          pd.Cevec[j+i*N] = pd.evec[j+i*N];
        }
      }
    }

    // invert rate matrix eigenvectors
    memcpy(pd.Cinv_evec,pd.Cevec,N2*sizeof(std::complex<double>));
    matrix_invert<std::complex<double> >(pd.Cinv_evec,N);
  }

  delete [] eval_cplx;
}

#endif


// get matrix exponential
//
void
matrix_exp_diag_scale(double* out_mat,  // dim*dim
                      const double scale,
                      matrix_exp_diag_prepdata& pd,
                      const bool is_get_deriv){

  const unsigned N = pd.dim;
  const unsigned N2 = N*N;

  if( ! pd.is_complex ){
    double* exp_eval = new double[N];

    for(unsigned i(0);i<N;++i){
      exp_eval[i] = std::exp(pd.eval[i]*scale);
      if(is_get_deriv) exp_eval[i] *= pd.eval[i];
    }

    // mult eigenvector by diag(exp(eval))
    double* evec_scale = new double[N2];
    for(unsigned i(0);i<N;++i){
      const double s=exp_eval[i];
      for(unsigned j(0);j<N;++j){
        evec_scale[j+i*N] = pd.evec[j+i*N]*s;
      }
    }

    delete [] exp_eval;

    // inv eigenvector * eigenvector matrix-matrix mult
    matrix_mult(pd.inv_evec,evec_scale,out_mat,N);

    delete [] evec_scale;

  } else {
    std::complex<double>* exp_eval = new std::complex<double>[N];

    for(unsigned i(0);i<N;++i){
      exp_eval[i] = std::exp(pd.Ceval[i]*scale);
      if(is_get_deriv) exp_eval[i] *= pd.Ceval[i];
    }

    // mult eigenvector by diag(exp(eval))
    std::complex<double>* evec_scale = new std::complex<double>[N2];
    for(unsigned i(0);i<N;++i){
      const std::complex<double>& s=exp_eval[i];
      for(unsigned j(0);j<N;++j){
        evec_scale[j+i*N] = pd.Cevec[j+i*N]*s;
      }
    }

    delete [] exp_eval;

    std::complex<double>* mat_exp = new std::complex<double>[N2];
    // inv eigenvector * eigenvector matrix-matrix mult
    matrix_mult(pd.Cinv_evec,evec_scale,mat_exp,N);

    delete [] evec_scale;

    // reformat matrix exponential for output
    for(unsigned i(0);i<N2;++i){
      out_mat[i] = std::real(mat_exp[i]);
      if(std::fabs(std::imag(mat_exp[i]))>MATRIX_EXP_ZERO_THRESH){
        log_os << "imaginary matrix_exp output: " << mat_exp[i] << "\n";
        log_os << "scale: " << scale << "\n";
        throw math_util_exception("imaginary matrix_exp output");
      }
    }

    delete [] mat_exp;
  }

  /// !!!!!!!!!!!!!non-portable sloppy hack!!!!!!!!!!!!!!!
  // assume we're writing out probabilities,so round very
  // small negative numbers to 0
  for(unsigned i(0);i<N2;++i){
    if(out_mat[i]<0. &&
       -out_mat[i]<MATRIX_EXP_NEG_THRESH) {
#ifdef DEBUG
      log_os << "WARNING:: matrix_diag_exp: zeroing:: " << out_mat[i] << "\n";
#endif
      out_mat[i] = 0.;
    }
  }
}

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

// $Id: matrix_invert.hh 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "clapack_util.h"
#include "math_util_exception.h"

#include <complex>

namespace CLAPACK {
  // overloaded c++ interfaces to lapack:
  inline
  void
  cpp_getrf(integer* m,
            integer* n,
            double* a,
            integer* lda,
            integer* pivots,
            integer* info){
    dgetrf_(m,n,a,lda,pivots,info);
  }

  inline
  void
  cpp_getrf(integer* m,
            integer* n,
            std::complex<double>* a,
            integer* lda,
            integer* pivots,
            integer* info){
    zgetrf_(m,n,reinterpret_cast<doublecomplex*>(a),lda,pivots,info);
  }

  inline
  void
  cpp_getri(integer* n,
            double* a,
            integer* lda,
            integer* pivots,
            double* work,
            integer* lwork,
            integer* info){
    dgetri_(n,a,lda,pivots,work,lwork,info);
  }

  inline
  void
  cpp_getri(integer* n,
            std::complex<double>* a,
            integer* lda,
            integer* pivots,
            std::complex<double>* work,
            integer* lwork,
            integer* info){
    zgetri_(n,reinterpret_cast<doublecomplex*>(a),lda,pivots,
            reinterpret_cast<doublecomplex*>(work),lwork,info);
  }

  inline
  double
  realsize(const double a) { return a; }

  inline
  double
  realsize(const std::complex<double> a) { return std::real(a); }
}


// invert matrix inplace
//
template <typename FloatType>
void
matrix_invert(FloatType* mat,
              const unsigned N){

  using namespace CLAPACK;

  integer n(static_cast<integer>(N));
  integer* pivots(new integer[N]);
  integer info(0);

  // LU decomposition
  cpp_getrf(&n,&n,mat,&n,pivots,&info);

  if(info){ throw math_util_exception("getrf (gen lu decomp) failure."); }

  FloatType dummy_work;
  integer work_size(-1);

  // workspace query
  cpp_getri(&n,mat,&n,pivots,&dummy_work,&work_size,&info);

  if(info){ throw math_util_exception("getri (matrix inverse) worksize failure in matrix_exp_diag_invert"); }

  work_size = static_cast<integer>(realsize(dummy_work));
  FloatType* work(new FloatType[work_size]);

  // inversion
  cpp_getri(&n,mat,&n,pivots,work,&work_size,&info);

  if(info){ throw math_util_exception("getri (matrix inverse) failure in matrix_exp_diag"); }

  delete[] pivots;
  delete[] work;
}

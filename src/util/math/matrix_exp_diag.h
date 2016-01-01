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

// $Id: matrix_exp_diag.h 743 2007-08-14 15:47:12Z ctsa $

/// \file
///
/// \brief matrix exponentiation by eigenvalue decomposition
///

#ifndef __MATRIX_EXP_DIAG_H
#define __MATRIX_EXP_DIAG_H


#include <complex>


struct matrix_exp_diag_prepdata;


/// \brief single step matrix exp by diag. not efficient for multi time calls.
///
template <typename FloatType>
void
matrix_exp_diag(FloatType* out_mat,       // dim*dim
                const FloatType* in_mat,  // dim*dim
                const FloatType scale,
                const unsigned N);

#ifdef USE_LAPACK
// get matrix exponential eigenvectors
//
void
matrix_exp_diag_prep(matrix_exp_diag_prepdata& pd,
                     const double* in_mat);  // dim*dim
void
matrix_exp_diag_prep(matrix_exp_diag_prepdata& pd,
                     const float* in_mat);   // dim*dim
#endif

// get matrix exponential
//
void
matrix_exp_diag_scale(double* out_mat,    // dim*dim
                      const double scale,
                      matrix_exp_diag_prepdata& pd,
                      const bool is_get_deriv=false);
void
matrix_exp_diag_scale(float* out_mat,     // dim*dim
                      const float scale,
                      matrix_exp_diag_prepdata& pd,
                      const bool is_get_deriv=false);


struct matrix_exp_diag_prepdata {

  matrix_exp_diag_prepdata(const unsigned _dim) : dim(_dim) {
    eval = new double[dim+2*(dim*dim)];
    evec = eval+dim;
    inv_evec = evec+dim*dim;
    Ceval = new std::complex<double>[dim+2*(dim*dim)];
    Cevec = Ceval+dim;
    Cinv_evec = Cevec+dim*dim;
  }

  ~matrix_exp_diag_prepdata(){
    delete [] eval;
    delete [] Ceval;
  }

  unsigned dim;
  bool is_complex;
  double* eval;
  double* evec;
  double* inv_evec;
  std::complex<double>* Ceval;
  std::complex<double>* Cevec;
  std::complex<double>* Cinv_evec;
};

#include "matrix_exp_diag.hh"

#endif

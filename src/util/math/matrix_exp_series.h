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

// $Id: matrix_exp_series.h 1055 2007-12-08 04:27:31Z ctsa $

/// \file
///
/// \brief matrix exponential by Taylor series approximation
///

#ifndef __MATRIX_EXP_SERIES_H
#define __MATRIX_EXP_SERIES_H

#include "math_util_exception.h"
#include "matrix_util.h"
#include "../general/workspace.h"


// series continues until all matrix values are stable (but not
// necessarily accurate) to this precision:
//
// this precision was chosen so that the final likelihood of a simple
// codon model using the swanson/yang bigmhc alignment matched the
// likelihood of the expm-diag model within < 1.25e-4 lnP
//
const double MATRIX_EXP_ABS_PREC(1e-12);

template <typename FloatType> struct matrix_exp_series_prepdata;

/// \brief single step matrix exp by series
///
template <typename FloatType>
void
matrix_exp_series(FloatType* out_mat,        // dim*dim
                  const FloatType* in_mat,   // dim*dim
                  const FloatType in_scale,
                  const unsigned dim);


/// \brief multi-step (part 1) matrix exp by series
///
template <typename FloatType>
void
matrix_exp_series_prep(matrix_exp_series_prepdata<FloatType>& mepd,
                       const FloatType* in_mat);    // dim*dim

/// \brief multi-step (part 2) matrix exp by series
///
// throws math_util_exception
template <typename FloatType>
void
matrix_exp_series_scale(FloatType* out_mat,          // dim*dim
                        matrix_exp_series_prepdata<FloatType>& mepd,
                        const FloatType in_scale,
                        workspace<char>& ws);


/// convenience struct for matrix_exp_series data:
///
template <typename FloatType>
struct matrix_exp_series_prepdata {

  struct iter_info {
    FloatType max;
    bool is_calc;
  };

  matrix_exp_series_prepdata(const unsigned _dim,
                             workspace<char>& ws1_init,
                             workspace<char>& ws2_init)
    : dim(_dim), size(0), prep(0), info(0), maxd(1.), is_symm(false),
      ws1(ws1_init), ws2(ws2_init) { grow(); }

  void grow();

  // data:
  unsigned dim;
  unsigned size;
  FloatType* prep;
  iter_info* info;
  FloatType maxd; // scaling factor
  bool is_symm;
  workspace<char>& ws1;
  workspace<char>& ws2;

  enum { CHUNK_SIZE = 25 };

  // set a hard series limit, well past the ridiculous point:
  enum { SERIES_MAX = 350 };
};


#include "matrix_exp_series.hh"

#endif

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

// $Id: matrix_exp_diag_float_thresh.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __MATRIX_EXP_DIAG_FLOAT_THRESH_H
#define __MATRIX_EXP_DIAG_FLOAT_THRESH_H


// hack for prob models... if neg numbers are less than this value,
// we push them to zero:
//
const double MATRIX_EXP_NEG_THRESH = 1.e-12;

// used to check that imaginary component is zero when expected:
const double MATRIX_EXP_ZERO_THRESH = 1.e-11;

// below this difference eigenvalues are considered equal
const double MATRIX_EXP_SAME_EVAL_THRESH = 1.e-8;

// below this difference return values are considered equal
const double MATRIX_EXP_SAME_RETURN_THRESH = 1.e-8;

#endif

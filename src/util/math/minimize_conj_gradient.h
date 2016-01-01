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

// $Id: minimize_conj_gradient.h 1055 2007-12-08 04:27:31Z ctsa $

/// \file

#ifndef __MINIMIZE_CONJ_GRADIENT_H
#define __MINIMIZE_CONJ_GRADIENT_H

#include "minfunc_interface.h"

/// \brief polak-ribiere conjugate gradient minimizer
///
/// \param vec_current starting function argument; on exit set to
/// minimum function argument
///
/// \param mf minimization function; must fulfill
/// minfunc_interface<FloatType>
///
/// \param tol final tolerance (on min function value improvement in
/// each iteration), minimization terminates when min function
/// improvement does not exceed this value in one iteration
///
/// \param f_min final minimum function value
///
/// \param iter number of iterations required to reach minimum
///
/// \param final_iter_delta_f decrease in f_min in the final iteration
///
/// \param max_iter maximum number of conj dir iterations to attempt
///
template <typename FloatType,
          typename MinFunc> // = minfunc_gradient_interface<FloatType> >
void
minimize_conj_gradient(FloatType* vec_current,
                       MinFunc& mf,
                       const FloatType tol,
                       FloatType& f_min,
                       unsigned& iter,
                       FloatType& final_iter_delta_f,
                       const unsigned max_iter);


#include "minimize_conj_gradient.hh"

#endif

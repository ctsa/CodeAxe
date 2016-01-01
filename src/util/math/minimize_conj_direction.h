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

// $Id: minimize_conj_direction.h 1055 2007-12-08 04:27:31Z ctsa $

/// \file

#ifndef __MINIMIZE_CONJ_DIRECTION_H
#define __MINIMIZE_CONJ_DIRECTION_H

#include "minfunc_interface.h"

/// \brief minimize function using modified conjugate direction set
/// method from (Powell, 1964)
///
/// \param vec_current starting function argument; on exit set to
/// minimum function argument
///
/// \param conj_direction_set starting conjugate direction set, scale
/// of conjugate vectors provides a hint to the minimizer of the scale
/// of parameter space to explore
///
/// \param mf minimization function; must fulfill
/// minfunc_interface<FloatType>
///
/// \param start_tol starting tolerance (on min function value
/// improvement in each iteration) this will be decreased by one order
/// of magnitude at each iteration where tolerance is not exceeded
/// until it reaches (or falls below) end_tol, at which point
/// minimization will terminate when end_tol is not exceeded. This is
/// a heuristic speed-up for long minimizations.
///
/// \param end_tol final tolerance (on min function value improvement
/// in each iteration)
///
/// \param line_tol minimum change in argument value per iteration during
/// line minimization
///
/// \param f_min final minimum function value
///
/// \param iter number of iterations required to reach minimum
///
/// \param final_iter_delta_f function improvement in final iteration
///
/// \param max_iter maximum number of conj dir iterations to attempt
///
template <typename FloatType,
          typename MinFunc> // = minfunc_interface<FloatType>
void
minimize_conj_direction(FloatType* vec_current,
                        FloatType* conj_direction_set,
                        MinFunc& mf,
                        FloatType& start_tol,
                        const FloatType end_tol,
                        const FloatType line_tol,
                        FloatType& f_min,
                        unsigned& iter,
                        FloatType& final_iter_delta_f,
                        const unsigned max_iter);

#include "minimize_conj_direction.hh"

#endif

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

// $Id: minimize_1d.h 880 2007-10-02 01:51:54Z ctsa $

/// \file

#ifndef __MINIMIZE_1D_H
#define __MINIMIZE_1D_H

//#include "minfunc_interface.h"


/// \brief Find local minimum in a 1d function
///
/// f1 must be set to mf.val(x1) at call
/// x2 is the first argument evaluated in the search for the minimum --
///    this is useful for suggesting the range of the 1d search
///
/// default minimization tolerance (x_tol) is sqrt(epsilon(FloatType))
/// as suggested by RP Brent in "Minimization without derivatives"
///
/// an optional secondary tolerance on the change to the function
/// value at each iteration may also be set
///
/// \param x1 starting function arg value
///
/// \param x2 next function arg to evaluate. abs(x1-x2) suggests the
/// scale of the domain to search over -- this can speed up the
/// minimization if known
///
/// \param f1 must be set to f(x1) [mf.val(x1)]
///
/// \param mf the minimization function. Must conforms to the
/// interface of: minfunc_1d_interface<FloatType> >
///
/// \param x_min returns the argument of minimum function value
///
/// \param f_min returns the minimum function value
///
/// \param x_tol "normal" minimization tolerance: the minimum change
/// in x per iteration
///
/// \param f_tol secondary minimization tolerance: the minimum change
/// in f(x) per iteration, if both tolerances are set the routine
/// stops when either criteria is not met in a single iteration
///
template <typename FloatType,
          typename MinFunc> // = minfunc_1d_interface<FloatType> >
void
minimize_1d(FloatType x1,
            FloatType x2,
            FloatType f1,
            MinFunc& mf,
            FloatType& x_min,
            FloatType& f_min,
            FloatType x_tol=0,
            FloatType f_tol=0);


#include "minimize_1d.hh"

#endif

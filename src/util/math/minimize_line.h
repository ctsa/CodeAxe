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

// $Id: minimize_line.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __MINIMIZE_LINE_H
#define __MINIMIZE_LINE_H


#include "minfunc_interface.h"

/// \brief minimize function along a line in a multidimensional parameter space
///
/// \param vec_base min line defined by vec_base+C*vec_dir; returns minimum
/// function argument
///
/// \param vec_dir min line defined by vec_base+C*vec_dir; length of
/// vec_dir gives minimizer a hint of the scale to search over
///
/// \param f_base must be set to f(vec_base) [ mf.val(vec_base ]
///
/// \param mf minimization function; must fulfill
/// minfunc_interface<FloatType>
///
/// \param f_min minimum function value
///
/// \param workspace tmp storage for routine: size=mf.dim()
///
/// \param x_tol minimum change in argument value per iteration
///
/// \param f_tol minimum change to function value per iteration
///
template <typename FloatType,
          typename MinFunc> // = minfunc_interface<FloatType> >
void
minimize_line(FloatType* vec_base,
              FloatType* vec_dir,
              FloatType f_base,
              MinFunc& mf,
              FloatType& f_min,
              FloatType* workspace,
              FloatType x_tol=0,
              FloatType f_tol=0);


#include "minimize_line.hh"

#endif

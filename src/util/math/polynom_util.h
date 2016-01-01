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

// $Id: polynom_util.h 1000 2007-11-05 05:06:18Z ctsa $

/// \file

#ifndef __POLYNOM_UTIL_H
#define __POLYNOM_UTIL_H

#include <algorithm>


/// \brief simple quadratic solver
///
/// parameters: a*x^2+b*x+c = 0
///
/// solution:    x1 >= x2
///
/// note no checks against (sqrt(-x) && FloatType != complex_type)
///
template <typename FloatType>
void
solve_quadratic(const FloatType a,
                const FloatType b,
                const FloatType c,
                FloatType& x1,
                FloatType& x2){

  const FloatType tmp1(-b/(2*a));
  const FloatType tmp2(std::sqrt(b*b-4*a*c)/(2*a));
  x1 = tmp1+tmp2;
  x2 = tmp1-tmp2;
  if(x1 < x2) std::swap(x1,x2);
}


#endif

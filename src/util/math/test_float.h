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

// $Id: test_float.h 743 2007-08-14 15:47:12Z ctsa $

/// \file
///
/// \brief machine portable floating point tests
///

#ifndef __TEST_FLOAT_H
#define __TEST_FLOAT_H

// these tests not working??? check against -ffast-math (and friends) in gcc
//
// add this if sassert is available: STATIC_ASSERT( (0./0.) != (0./0.) )


// replacements for BSD libc's
#ifdef DARWIN_HACK

template <typename FloatType>
bool
is_float_nan(FloatType x){ return x != x; }

template <typename FloatType>
bool
is_float_inf(FloatType x){ return is_float_nan(x*0) && (! is_float_nan(x)); }

#else

#include <cmath>

template <typename FloatType>
bool
is_float_nan(FloatType x){ return isnan(x); }

template <typename FloatType>
bool
is_float_inf(FloatType x){ return isinf(x); }

#endif


template <typename FloatType>
bool
is_float_invalid(FloatType x){
  return is_float_nan(x) || is_float_inf(x);
}

#endif


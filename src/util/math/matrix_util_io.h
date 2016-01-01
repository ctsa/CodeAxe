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

// $Id: matrix_util_io.h 743 2007-08-14 15:47:12Z ctsa $

/// \file
///
/// \brief simple matrix i/o
///

/// tmp hack to keep irv_t<T> overload working correctly
/// ...locking for a real solution
#include "indy_random_var_io.h"


#ifndef __MATRIX_UTIL_IO_H
#define __MATRIX_UTIL_IO_H


#include <iosfwd>

enum mutil_bogus_t { MAX_MATRIX_PRINT = 21 };


/// \brief generic NxN matrix printer
///
template <typename FloatType>
void
matrix_report(const FloatType mat[],
              const unsigned N,
              const char* syms,
              std::ostream& os,
              const unsigned cell_width = 8,
              const unsigned prec = 2);



/// \brief generic NxN matrix debug printer
///
template <typename FloatType>
void
matrix_dump(const FloatType mat[],
            const unsigned N,
            std::ostream& os,
            int prec = 2);


#include "matrix_util_io.hh"

#endif

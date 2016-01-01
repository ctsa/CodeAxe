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

// $Id: linear_solver.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __LINEAR_SOLVER_H
#define __LINEAR_SOLVER_H

#ifdef USE_LAPACK

/// \brief solve the system of linear equations Ax=b for x by LU
/// decomposition of A
///
void
linear_solver(const double* A,
              const double* b,
              double* x,
              const unsigned N);

/// \brief solve the system of linear equations Ax=b for x by LU
/// decomposition of A. x is improved by iterative refinement
///
/// \param ferr upper bound on the largest error term in x
///
void
linear_solver_refined(const double* A,
                      const double* b,
                      double* x,
                      double& ferr,
                      const unsigned N);

#endif

#endif

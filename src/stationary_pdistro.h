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

// $Id: stationary_pdistro.h 744 2007-08-14 18:09:31Z ctsa $

/// \file
///
/// \brief utilities for calculating the stationary distribution from
/// rate/trans prob matrices
///

#ifndef __STATIONARY_PDISTRO_H
#define __STATIONARY_PDISTRO_H

#include "subs_ml_types.h"

/// \todo generalize non-iterative solution and move into util/math

/// \brief get the stationary distribution of the rate matrix
///
/// gets the stationary distribution by LU decomposition. if lapack is
/// not available, falls back to an iterative solution, which starts
/// from a uniform distribution unless flag is set to start from
/// stat_pdistro's input value
///
///  throws substk_exception
void
get_stationary_pdistro_from_rates(smlfloat* stat_pdistro,
                                  const smlfloat* rates,
                                  const unsigned N,
                                  const bool is_start_from_input = false);

/// throws substk_exception
void
get_stationary_pdistro_from_pmatrix(smlfloat* stat_pdistro,
                                    const smlfloat* pstep,
                                    const unsigned N,
                                    const bool is_start_from_input = false);

#endif

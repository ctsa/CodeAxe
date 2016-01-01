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

// $Id: sim_discrete.h 1073 2008-01-03 01:35:05Z ctsa $

/// \file

#ifndef __SIM_DISCRETE_H
#define __SIM_DISCRETE_H

#include "nuc_seq_data.h"
#include "sim_options.h"
#include "subs_ml_model.h"

/// \brief non-context discrete time sim
///
/// this function has a subset of the capabilities of the context-dependent
/// discrete simulator, so there's no reason for it to remain except as
/// a confidence check on the more general sim functions
///
///
void
simulate_data_discrete(const sim_options& sim_opt,
                       const subs_ml_model& mdl,
                       nuc_seq_data& nsd);

#endif

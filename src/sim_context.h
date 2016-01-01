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

// $Id: sim_context.h 1071 2008-01-02 22:58:00Z ctsa $

/// \file

#ifndef __SIM_CONTEXT_H
#define __SIM_CONTEXT_H

#include "nuc_seq_data.h"
#include "sim_options.h"
#include "subs_ml_model.h"


/// \brief continuous and discrete time data simulation for models
/// with context dependent mutation rates
///
void
simulate_data_context(const sim_options& sim_opt,
                      const subs_ml_model& mdl,
                      nuc_seq_data& nsd);

#endif

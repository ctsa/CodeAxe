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

// $Id: subs_prob_util.h 1107 2008-01-25 02:40:20Z ctsa $

/// \file

#ifndef __SUBS_PROB_UTIL_H
#define __SUBS_PROB_UTIL_H

#include "rates_func_options.h"
#include "subs_ml_types.h"
#include "util/general/workspace.h"

#include <vector>

struct subs_ml_model;

/// \brief get the transition prob matrices for every branch in the model
///
/// this is the malloc on first call only version.
///
void
get_all_branch_subs_prob(prob_t** branch_prob,
                         const subs_ml_model& mdl,
                         const rates_func_options_base& opt,
                         workspace<char> ws[3]);

/// \brief get the transition prob matrices for every branch in the model
///
inline
void
get_all_branch_subs_prob(prob_t** branch_prob,
                         const subs_ml_model& mdl,
                         const rates_func_options_base& opt){
  workspace<char> ws[3];
  get_all_branch_subs_prob(branch_prob,mdl,opt,ws);
}

/// \brief get a single transition prob matrix for the model
///
void
get_single_branch_subs_prob(prob_t* branch_prob,
                            const smlfloat time,
                            const subs_ml_model& mdl,
                            const rates_func_options_base& opt,
                            const unsigned branch_id,
                            workspace<char> ws[3]);

#endif

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

// $Id: subs_ml_model_root_from_obs.h 801 2007-09-08 19:13:25Z ctsa $

/// \file
///
/// \brief "parameterless" root distros for a model from observed leaf distros
///

#ifndef __SUBS_ML_MODEL_ROOT_FROM_OBS_H
#define __SUBS_ML_MODEL_ROOT_FROM_OBS_H

#include "subs_ml_model.h"

/// \brief time weighted average of observed state distros at root
///
void
subs_ml_model_obs_state_time_avg_root(const subs_ml_model& mdl,
                                      prob_t* root_down_prob,
                                      const unsigned cat_no);

#ifdef BRANCH_SPECIFC_OBS
bool
subs_ml_model_root_node_prob_nspround(const subs_ml_model& mdl,
                                      prob_t* root_node_prob,
                                      const unsigned cat_no,
                                      prob_t** node_state_prob);
#endif

/// \brief approximate distro of states at tree root using the
/// per-branch least squares ('BLS') method
///
/// root distro estimated from observed state distros at each leaf and
/// model rate matrix
///
bool
subs_ml_model_root_node_prob(const subs_ml_model& mdl,
                             prob_t* root_node_prob,
                             const unsigned cat_no);

#endif

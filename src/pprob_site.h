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

// $Id: pprob_site.h 1161 2008-03-20 19:28:26Z ctsa $

/// \file

#ifndef __PPROB_SITE_H
#define __PPROB_SITE_H

#include "bi_tree.h"
#include "site_code.h"
#include "subs_ml_types.h"
#include "tree_probs.h"

/// \brief get interior node posterior probs for site
///
/// also returns leaf node posterior as a trivial result, it could be
/// meaningful at some point with trickier missing/gap state handlers.
///
void
get_site_pprob_single_cat(const site_code::index_type* seq_state,
                          const bi_tree& t,
                          const tree_probs& tprob,
                          const unsigned n_states,
                          prob_t** pprob); // n_nodes*n_states



/// \brief get interior node posterior probs for site (no-malloc version)
///
void
get_site_pprob_single_cat_prep(const site_code::index_type* seq_state,
                               const bi_tree& t,
                               const tree_probs& tprob,
                               const unsigned n_states,
                               prob_t** up_branch_prob, // n_branches*n_states
                               prob_t** down_branch_prob, // n_branches*n_states
                               prob_t* tmp_state_array,  // n_states
                               prob_t** pprob); // n_nodes*n_states


/// \brief get interior node posterior probs for site category mixture model
///
void
get_site_pprob_from_branch_prob(const bi_tree& tree,
                                const unsigned n_states,
                                const prob_t* root_prob,
                                const prob_t* const * up_branch_prob,
                                const prob_t* const * down_branch_prob,
                                const site_code::index_type* seq_state,
                                prob_t** pprob);

#endif

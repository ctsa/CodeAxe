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

// $Id: lhood_site.h 1161 2008-03-20 19:28:26Z ctsa $

/// \file
/// \brief pruning algorithm site probability functions
///

#ifndef __LHOOD_SITE_H
#define __LHOOD_SITE_H

#include "bi_tree.h"
#include "site_code.h"
#include "subs_ml_types.h"
#include "tree_probs.h"

/// !!ALL TPROBS TRANSPOSED IN ALL FUNCTIONS!!


/// \brief site prob for one category by pruning
///
prob_t
get_site_prob_single_cat(const site_code::index_type* seq_state,
                         const bi_tree& t,
                         const tree_probs& tprob,
                         const unsigned n_states);

/// \brief site prob for one category by pruning (no-malloc version)
///
/// note that down_prob_tmp size is suboptimal, but isn't worth
/// worrying about for small trees
///
prob_t
get_site_prob_single_cat_prep(const site_code::index_type* seq_state,
                              const bi_tree& t,
                              const tree_probs& tprob,
                              const unsigned n_states,
                              prob_t* const * const down_prob_tmp); // n_branches*n_states

/// this provides, for each value of root_pp[i] the probability of the
/// observed site if the root node is i, and the prior probability of
/// state i at the root is 1
///
///
void
get_site_root_partial_prob_single_cat_prep(const site_code::index_type* seq_state,
                                           const bi_tree& t,
                                           const tree_probs& tprob,
                                           const unsigned n_states,
                                           prob_t* const root_pp,
                                           prob_t* const * const down_prob_tmp); // n_branches*n_states

/// \brief category mixture model site prob by pruning
///
prob_t
get_site_prob(const site_code::index_type* seq_state,
              const bi_tree& t,
              const tree_probs* tprob,
              const prob_t* site_cat_prior_pdistro,
              const unsigned n_states,
              const unsigned n_site_cats);

#endif

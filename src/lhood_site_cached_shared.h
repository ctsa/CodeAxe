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

// $Id: lhood_site_cached_shared.h 743 2007-08-14 15:47:12Z ctsa $

/// \file
/// \brief components shared between regular and em lhood_site_cached calculations
///

#ifndef LHOOD_SITE_CACHED_SHARED_H
#define LHOOD_SITE_CACHED_SHARED_H

#include "lhood_site_cached_tree_lhood_info.h"
#include "simple_util.h"
#include "subs_ml_types.h"


struct tree_lhood_cache {
  tree_lhood_cache(const bi_tree_node* n=0,
                   const prob_t* v=0)
    : node(n), val(v) {}

  const bi_tree_node* node;
  const prob_t* val;
};


/// \brief stores intermediate results from partial tree evaluation
///
/// branch_lhood.val: [n_states] partial results for tree branch extending from
///                   branch_lhood.node to the unevaluated node
///
struct tree_walk_cache {

  tree_walk_cache(unsigned n)
    : branch_lhood_alloced(new prob_t[n]),
      combined_lhood_alloced(new prob_t[n]),
      branch_lhood(0,branch_lhood_alloced),
      combined_lhood(0,combined_lhood_alloced){}

  ~tree_walk_cache(){
    if(branch_lhood_alloced) delete [] branch_lhood_alloced;
    if(combined_lhood_alloced) delete [] combined_lhood_alloced;
  }

  prob_t* branch_lhood_alloced;
  prob_t* combined_lhood_alloced;

  tree_lhood_cache branch_lhood;
  tree_lhood_cache combined_lhood;
};



/// info gathered from all org levels. for each site record the cache
/// for the branch above each node in both the up and down direction,
/// for all nodes where this info is required for the em estep:
///
/// this structure is used both for EM and internal node post prob
/// calculations
///
struct tree_walk_cache_estep_top {

  tree_walk_cache_estep_top(const unsigned b,
                            const unsigned s)
    : up_branch_prob(b,s),down_branch_prob(b,s),
      site_prob(0.) {}

  simple_matrix<prob_t> up_branch_prob;
  simple_matrix<prob_t> down_branch_prob;
  prob_t site_prob;
};



/// \brief entry point for each tree leaf into the tree lhood
/// calculation
///
/// evaluation will proceed as far as possible from this leaf, given
/// the information from previously evaluated leaves supplied in
/// prior_branch_lhood, the partial lhood calculation is returned in
/// tw
///
void
lhood_site_cached_combine_branches(const tree_lhood_info& ti,
                                   const bi_tree_node* current_node,
                                   std::vector<tree_lhood_cache>& prior_branch_lhood,
                                   tree_walk_cache& tw,
                                   tree_walk_cache_estep_top* twe = 0);

#endif

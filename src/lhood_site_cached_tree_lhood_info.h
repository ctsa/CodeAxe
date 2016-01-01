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

// $Id: lhood_site_cached_tree_lhood_info.h 916 2007-10-12 21:13:20Z ctsa $

/// \file

#ifndef __TREE_LHOOD_INFO_H
#define __TREE_LHOOD_INFO_H

#include "bi_tree.h"
#include "tree_probs.h"

/// \brief all const info used in lhood calc
///
struct tree_lhood_info {
  tree_lhood_info(const bi_tree& init_tree,
                  const tree_probs& init_tprob,
                  const int init_n_states,
                  const int init_ambiguous_leaf_state = -1)
    : tree(init_tree), tprob(init_tprob), n_states(init_n_states) {
    if(init_ambiguous_leaf_state==-1) ambiguous_leaf_state = n_states;
    else                              ambiguous_leaf_state = init_ambiguous_leaf_state;
  }


  const bi_tree& tree;
  const tree_probs& tprob;
  const int n_states;
  unsigned ambiguous_leaf_state;
};


#endif

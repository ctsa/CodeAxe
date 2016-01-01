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

// $Id: lhood_site.cc 1161 2008-03-20 19:28:26Z ctsa $

/// \file

#include "lhood_site.h"
#include "util/math/matrix_util.h"
#include "util/general/die.h"


// probability of subtree rooted at node b, for all n_states values
// at node b:
//
static
void
tree_branch_down_prob(const bi_tree_node* b,
                      const tree_probs& tprob,
                      const site_code::index_type* seq_state,
                      const unsigned n_states,
                      const prob_t*& branch_p,
                      prob_t* const * const down_prob_tmp){

  const unsigned branch_id(b->branch_id());

  if(b==0){
    die("tree_branch_down_prob: null tree node");

  } else if(b->is_leaf()){  // leaf node
    const unsigned seq_val(seq_state[b->leaf_id()]);
    branch_p = tprob.branch_prob(branch_id)+seq_val*n_states;

  } else {  // assume child2 exists
    const prob_t* child1_branch_p(0);
    const prob_t* child2_branch_p(0);

    tree_branch_down_prob(b->child1(),tprob,seq_state,n_states,
                          child1_branch_p,down_prob_tmp);
    tree_branch_down_prob(b->child2(),tprob,seq_state,n_states,
                          child2_branch_p,down_prob_tmp);

    prob_t* branch_p_tmp(down_prob_tmp[branch_id]);
    for(unsigned i_parent(0);i_parent<n_states;++i_parent){ branch_p_tmp[i_parent]=0.; }

    for(unsigned i_node(0);i_node<n_states;++i_node){
      const prob_t node_down_p(child1_branch_p[i_node]*child2_branch_p[i_node]);
      const prob_t* const tprob_vec(tprob.branch_prob(branch_id)+i_node*n_states);
      for(unsigned i_parent(0);i_parent<n_states;++i_parent){
        branch_p_tmp[i_parent] += tprob_vec[i_parent]*node_down_p;
      }
    }

    branch_p = branch_p_tmp;

  }
}



void
get_site_root_partial_prob_single_cat_prep(const site_code::index_type* seq_state,
                                           const bi_tree& t,
                                           const tree_probs& tprob,
                                           const unsigned n_states,
                                           prob_t* const root_pp,
                                           prob_t* const * const down_prob_tmp){ // n_branches*n_states

  if(t.root()->is_leaf()) {
    const prob_t up(1./static_cast<prob_t>(n_states));
    for(unsigned i(0);i<n_states;++i) root_pp[i] = up;
    return;
  }

  const prob_t* child1_branch_p(0);
  const prob_t* child2_branch_p(0);

  tree_branch_down_prob(t.root()->child1(),tprob,seq_state,n_states,
                        child1_branch_p,down_prob_tmp);
  tree_branch_down_prob(t.root()->child2(),tprob,seq_state,n_states,
                        child2_branch_p,down_prob_tmp);

  for(unsigned i(0);i<n_states;++i){
    root_pp[i] = child1_branch_p[i]*child2_branch_p[i];
  }
}



prob_t
get_site_prob_single_cat_prep(const site_code::index_type* seq_state,
                              const bi_tree& t,
                              const tree_probs& tprob,
                              const unsigned n_states,
                              prob_t* const * const down_prob_tmp){ // n_branches*n_states

  if(t.root()->is_leaf()) return 1.;

  const prob_t* child1_branch_p(0);
  const prob_t* child2_branch_p(0);

  tree_branch_down_prob(t.root()->child1(),tprob,seq_state,n_states,
                        child1_branch_p,down_prob_tmp);
  tree_branch_down_prob(t.root()->child2(),tprob,seq_state,n_states,
                        child2_branch_p,down_prob_tmp);

  prob_t p(0.);
  for(unsigned i(0);i<n_states;++i){
    p += tprob.root_prob()[i]*child1_branch_p[i]*child2_branch_p[i];
  }
  return p;
}




prob_t
get_site_prob_single_cat(const site_code::index_type* seq_state,
                         const bi_tree& t,
                         const tree_probs& tprob,
                         const unsigned n_states){

  simple_matrix<prob_t> down_prob_tmp(t.branch_size(),n_states);

  return get_site_prob_single_cat_prep(seq_state,t,tprob,n_states,down_prob_tmp.ptr());
}




prob_t
get_site_prob(const site_code::index_type* seq_state,
              const bi_tree& t,
              const tree_probs* tprob,
              const prob_t* site_cat_prior_pdistro,
              const unsigned n_states,
              const unsigned n_site_cats){

  prob_t p(0.);
  for(unsigned sc(0);sc<n_site_cats;++sc){
    p += site_cat_prior_pdistro[sc]*get_site_prob_single_cat(seq_state,t,
                                                             tprob[sc],n_states);
  }

  return p;
}

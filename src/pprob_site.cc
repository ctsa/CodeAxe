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

// $Id: pprob_site.cc 1161 2008-03-20 19:28:26Z ctsa $

/// \file

#include "pprob_site.h"
#include "substk_exception.h"
#include "util/math/matrix_util.h"
#include "util/math/prob_util.h"

#ifdef XTRA_DEBUG
#include "lhood_site.h"
#include "util/general/log.h"
#include <ostream>
#endif



static
void
tree_branch_down_prob(const bi_tree_node* b,
                      const tree_probs& tprob,
                      const unsigned n_states,
                      const prob_t* const * up_branch_prob,
                      const prob_t* root_prob,
                      prob_t** down_branch_prob,
                      prob_t* const tmp_state_array){

  if(b==0 || b->is_root()) throw substk_exception("tree_branch_down_prob: unexpected tree");

  const int branch_id(b->branch_id());
  prob_t* branch_p(down_branch_prob[branch_id]);

  const bi_tree_node* parent(b->parent());
  const bi_tree_node* sister(b->sister());

  const int parent_branch_id(parent->branch_id());
  const int sister_branch_id(sister->branch_id());

  const prob_t* parent_p;
  if(parent->is_root()){ // use root as parent down branch
    parent_p=root_prob;
  } else {
    parent_p=down_branch_prob[parent_branch_id];
  }

  for(unsigned i_node(0);i_node<n_states;++i_node){
    const prob_t* const tprob_vec(tprob.branch_prob(branch_id)+i_node*n_states);

    prob_t* parent_node_combined_p(tmp_state_array);

    for(unsigned i_parent(0);i_parent<n_states;++i_parent){
      parent_node_combined_p[i_parent] = up_branch_prob[sister_branch_id][i_parent]*parent_p[i_parent];
    }

    // dot product sums over parent node states...
    branch_p[i_node] = array_dot(tprob_vec,parent_node_combined_p,n_states);
  }

  if(! b->is_leaf()){
    tree_branch_down_prob(b->child1(),tprob,n_states,up_branch_prob,
                          root_prob,down_branch_prob,tmp_state_array);
    tree_branch_down_prob(b->child2(),tprob,n_states,up_branch_prob,
                          root_prob,down_branch_prob,tmp_state_array);
  }
}




// probability of subtree rooted at node b, for all NSTATE values
// at node b:
//
static
void
tree_branch_up_prob(const bi_tree_node* b,
                    const tree_probs& tprob,
                    const site_code::index_type* seq_state,
                    const unsigned n_states,
                    prob_t** up_branch_prob){

  const int branch_id(b->branch_id());
  prob_t* branch_p(up_branch_prob[branch_id]);

  if(b==0 || b->is_root()) throw substk_exception("tree_branch_up_prob: unexpected_tree");

  if(b->is_leaf()){
    const unsigned seq_val(seq_state[b->leaf_id()]);
    for(unsigned i(0);i<n_states;++i){
      branch_p[i] = tprob.branch_prob(branch_id)[seq_val*n_states+i];
    }

  } else {
    const int child1_branch_id(b->child1()->branch_id());
    const int child2_branch_id(b->child2()->branch_id());

    tree_branch_up_prob(b->child1(),tprob,seq_state,n_states,up_branch_prob);
    tree_branch_up_prob(b->child2(),tprob,seq_state,n_states,up_branch_prob);

    std::fill(branch_p,branch_p+n_states,0.);

    // loop over node state:
    for(unsigned i_node(0);i_node<n_states;++i_node){
      const prob_t node_down_combined_p(up_branch_prob[child1_branch_id][i_node]*
                                        up_branch_prob[child2_branch_id][i_node]);
      const prob_t* const tprob_vec(tprob.branch_prob(branch_id)+i_node*n_states);
      // loop over state of parent node:
      for(unsigned i_parent(0);i_parent<n_states;++i_parent){
        branch_p[i_parent] += tprob_vec[i_parent]*node_down_combined_p;
      }
    }
  }
}




void
get_site_pprob_from_branch_prob(const bi_tree& tree,
                                const unsigned n_states,
                                const prob_t* root_prob,
                                const prob_t* const * up_branch_prob,
                                const prob_t* const * down_branch_prob,
                                const site_code::index_type* seq_state,
                                prob_t** pprob){

#ifdef DEBUG
    const unsigned n_branches(tree.branch_size());
    for(unsigned j(0);j<n_branches;++j){
      for(unsigned i(0);i<n_states;++i){
        assert(up_branch_prob[j][i] >= 0.);
        assert(down_branch_prob[j][i] >= 0.);
      }
    }
#endif

  const unsigned n_nodes(tree.node_size());

  for(unsigned j(0);j<n_nodes;++j){
    const bi_tree_node* jn(tree.node(j));
    prob_t* pn(pprob[j]);

    if(jn->is_root()){ // root node:
      for(unsigned i(0);i<n_states;++i){
        pn[i] =
          root_prob[i]*
          up_branch_prob[jn->child1()->branch_id()][i]*
          up_branch_prob[jn->child2()->branch_id()][i];
      }
#ifdef XTRA_DEBUG
      prob_t sum(0.);
      for(unsigned i(0);i<n_states;++i) sum += pn[i];
      log_os << "node pp(r): " << j << " " << sum << "\n";
#endif
    } else if(jn->is_leaf()){
      std::fill(pn,pn+n_states,0.);
      pn[seq_state[jn->leaf_id()]] = 1.;

    } else {
      for(unsigned i(0);i<n_states;++i){
        pn[i] =
          down_branch_prob[jn->branch_id()][i]*
          up_branch_prob[jn->child1()->branch_id()][i]*
          up_branch_prob[jn->child2()->branch_id()][i];
      }
#ifdef XTRA_DEBUG
      prob_t sum(0.);
      for(unsigned i(0);i<n_states;++i) sum += pn[i];
      log_os << "node ppsum: " << j << " " << sum << "\n";
#endif
    }

    pdistro_norm(pn,pn+n_states);
  }
}




void
get_site_pprob_single_cat_prep(const site_code::index_type* seq_state,
                               const bi_tree& t,
                               const tree_probs& tprob,
                               const unsigned n_states,
                               prob_t** up_branch_prob, // n_branches*n_states
                               prob_t** down_branch_prob, // n_branches*n_states
                               prob_t* tmp_state_array,  // n_states
                               prob_t** pprob){ // n_nodes*n_states

  // temporary structures should match what's found in twe in the cached version of this routine:
  // basically just need an up and down prob distro for every branch

#ifdef DEBUG
  // fill with bogus vals to check against branch skip...
  const unsigned n_branches(t.branch_size());
  const unsigned vi(n_branches*n_states);
  std::fill(up_branch_prob[0],up_branch_prob[0]+vi,-1);
  std::fill(up_branch_prob[0],up_branch_prob[0]+vi,-1);
#endif

  tree_branch_up_prob(t.root()->child1(),tprob,seq_state,n_states,up_branch_prob);
  tree_branch_up_prob(t.root()->child2(),tprob,seq_state,n_states,up_branch_prob);

  tree_branch_down_prob(t.root()->child1(),tprob,n_states,up_branch_prob,
                        tprob.root_prob(),down_branch_prob,tmp_state_array);
  tree_branch_down_prob(t.root()->child2(),tprob,n_states,up_branch_prob,
                        tprob.root_prob(),down_branch_prob,tmp_state_array);

  // share pprob calculation step from cached_site version....
  get_site_pprob_from_branch_prob(t,n_states,tprob.root_prob(),
                                  up_branch_prob,down_branch_prob,
                                  seq_state,pprob);
}




void
get_site_pprob_single_cat(const site_code::index_type* seq_state,
                          const bi_tree& t,
                          const tree_probs& tprob,
                          const unsigned n_states,
                          prob_t** pprob){ // n_nodes*n_states

  // temporary structures should match what's found in twe in the cached version of this routine:
  const unsigned n_branches(t.branch_size());
  simple_matrix<prob_t> up_branch_prob(n_branches,n_states);
  simple_matrix<prob_t> down_branch_prob(n_branches,n_states);
  simple_array<prob_t> tmp_state_array(n_states);

  get_site_pprob_single_cat_prep(seq_state,t,tprob,n_states,
                                 up_branch_prob.ptr(),down_branch_prob.ptr(),tmp_state_array.ptr(),
                                 pprob);
}

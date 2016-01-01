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

// $Id: lhood_site_cached_shared.cc 999 2007-11-02 18:48:00Z ctsa $

/// \file

#include "lhood_site_cached_shared.h"
#include "util/math/matrix_util.h"

#ifdef LSC_DEBUG
#include "util/general/log.h"

#include <ostream>
#endif

// evaluate the branch between current node and combined_lhood_node
static
void
eval_branch_lhood(const tree_lhood_info& ti,
                  const bi_tree_node* current_node,
                  std::vector<tree_lhood_cache>& prior_branch_lhood,
                  tree_walk_cache& tw,
                  tree_walk_cache_estep_top* twe){

  tw.branch_lhood.val=tw.branch_lhood_alloced;
  tw.branch_lhood.node=tw.combined_lhood.node;

  // down branch
  if(current_node->parent()==tw.combined_lhood.node) {
#ifdef LSC_DEBUG
    log_os << "INTREE: down branch\n";
#endif
    const unsigned branch_id(current_node->branch_id());
    if(ti.tprob.is_branch_tprob_identity(branch_id)){
      for(int i(0);i<ti.n_states;++i){
        tw.branch_lhood_alloced[i] = tw.combined_lhood.val[i];
      }
    } else {
      const prob_t* const branch_tprob(ti.tprob.branch_prob(branch_id));
      for(int i(0);i<ti.n_states;++i){
        const prob_t* const tprob_vec(branch_tprob+i*ti.n_states);
        tw.branch_lhood_alloced[i] = array_dot(tprob_vec,tw.combined_lhood.val,ti.n_states);
      }
    }

    if(twe){
      for(int i(0);i<ti.n_states;++i){
        twe->down_branch_prob[branch_id][i] = tw.branch_lhood.val[i];
      }
#ifdef LSC_DEBUG
      log_os << "INTREE: twe down branchid/node " << branch_id << " " << current_node << "\n";
#endif
    }
  }

  // up branch:
  else {
#ifdef LSC_DEBUG
    log_os << "INTREE: up branch: " << tw.combined_lhood.node << "\n";
#endif
    const unsigned branch_id(tw.combined_lhood.node->branch_id());
    if(ti.tprob.is_branch_tprob_identity(branch_id)){
      for(int i(0);i<ti.n_states;++i){
        tw.branch_lhood_alloced[i] = tw.combined_lhood.val[i];
      }
    } else {
      const prob_t* const branch_tprob(ti.tprob.branch_prob(branch_id));
      for(int i(0);i<ti.n_states;++i){ tw.branch_lhood_alloced[i] = 0.; }

      for(int i(0);i<ti.n_states;++i){
        const prob_t& node_val(tw.combined_lhood.val[i]);
        const prob_t* const tprob_vec(branch_tprob+i*ti.n_states);
        for(int j(0);j<ti.n_states;++j){
          tw.branch_lhood_alloced[j] += tprob_vec[j]*node_val;
        }
      }
    }

    if(twe){
      for(int i(0);i<ti.n_states;++i){
        twe->up_branch_prob[branch_id][i] = tw.branch_lhood.val[i];
      }
#ifdef LSC_DEBUG
      log_os << "INTREE: twe up branchid/node " << branch_id << " " << current_node << "\n";
#endif
    }
  }

#ifdef DEBUG
  for(unsigned i(0);i<static_cast<unsigned>(ti.n_states);++i) assert(tw.branch_lhood.val[i]>=0.);
#endif

  lhood_site_cached_combine_branches(ti,current_node,prior_branch_lhood,tw,twe);
}



// combine known branch information at current_node, if possible.
void
lhood_site_cached_combine_branches(const tree_lhood_info& ti,
                                   const bi_tree_node* current_node,
                                   std::vector<tree_lhood_cache>& prior_branch_lhood,
                                   tree_walk_cache& tw,
                                   tree_walk_cache_estep_top* twe){

  const prob_t* second_branch_lhood(0);
  const bi_tree_node* next_node(0);

  // root case:
  if(current_node->is_root()){
    second_branch_lhood=ti.tprob.root_prob();
    next_node=tw.branch_lhood.node->sister();

#ifdef LSC_DEBUG
    log_os << "INTREE:: nodefunc root case. next_node= " << next_node << "\n";
#endif
  // test for prior case:
  } else {

    const bi_tree_node* edge_node[3] = { current_node->parent(),
                                         current_node->child1(),
                                         current_node->child2() };

    // put 2 continuing edges first in edge_node list
    for(unsigned i(0);i<2;++i){
      if(edge_node[i] == tw.branch_lhood.node) {
        std::swap(edge_node[i],edge_node[2]);
      }
    }

    bool is_prior(false);
    for(unsigned i(0);i<prior_branch_lhood.size();++i){
      for(unsigned j(0);j<2;++j) {
        if(edge_node[j] == prior_branch_lhood[i].node){
          is_prior=true;
          second_branch_lhood=prior_branch_lhood[i].val;
          next_node=edge_node[(j+1)%2];
#ifdef LSC_DEBUG
          log_os << "INTREE:: nodefunc prior case. next_node= " << next_node << "\n";
#endif
          break;
        }
      }
      if(is_prior) {
        std::swap(prior_branch_lhood[i],prior_branch_lhood.back());
        prior_branch_lhood.pop_back();
        break;
      }
    }

    // null case:
    if(! is_prior) return;
  }

  tw.combined_lhood.node=current_node;

  for(int i(0);i<ti.n_states;++i){
    tw.combined_lhood_alloced[i] = second_branch_lhood[i]*tw.branch_lhood.val[i];
  }
#ifdef LSC_DEBUG
  log_os << "INTREE current_node: " << current_node << "\n";
  for(unsigned i(0);i<4;++i){
    log_os << "combined lhood: " << i << " " << tw.combined_lhood.val[i] << "\n";
  }
#endif

  // check that we're not moving into a leaf node -- else time to
  // return node_cache back down the branch stack and move to next
  // leaf:
  if(! next_node->is_leaf()) {
    eval_branch_lhood(ti,next_node,prior_branch_lhood,tw,twe);
  }
}

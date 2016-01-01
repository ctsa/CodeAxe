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

// $Id: lhood_site_cached.cc 962 2007-10-25 22:18:44Z ctsa $

/// \file

#include "lhood_site_cached.h"
#include "lhood_site_cached_shared.h"
#include "util/general/die.h"
#include "util/math/matrix_util.h"

#ifdef LSC_DEBUG
#include "lhood_site.h"
#include "util/general/log.h"

#include <ostream>
#endif




static
void
eval_last_leaf(const tree_lhood_info& ti,
               const site_data_fastlup_core& sdf,
               const tree_lhood_cache& prior_combined_lhood,
               const int leaf_index,
               site_data_fastlup_core::index_type* const leaf_state,
               int& site_id,
               const std::vector<bool>& site_prob_mask,
               prob_t* const site_prob){

  const unsigned leaf_id(sdf.order[leaf_index]);
  const bi_tree_node* leaf_node(ti.tree.leaf_node(leaf_id));
  const int branch_id(leaf_node->branch_id());
  const prob_t* const branch_tprob(ti.tprob.branch_prob(branch_id));

  assert(leaf_index+1 == static_cast<int>(ti.tree.leaf_size()));
  assert(leaf_node->parent()==prior_combined_lhood.node);

#ifdef LSC_DEBUG
  log_os << "last_leaf: leaf_node: " << leaf_node << "\n";
  log_os << "last_leaf: parent_node: " << leaf_node->parent << "\n";
#endif

  while(true) {
    leaf_state[leaf_id] = sdf.data_core[site_id].index[leaf_id];
#ifdef LSC_DEBUG
    log_os << "eval site_id: " << site_id << "\n";
    for(unsigned i(0); i<ti.tree.n_leaves(); ++i){
      log_os << "os: " << i << " " << ti.tree.leaf[i]->label << " " << leaf_state[i] << "\n";
    }
#endif

    if(leaf_state[leaf_id] == ti.ambiguous_leaf_state){ // ambiguous state case
      site_prob[site_id]=array_sum(prior_combined_lhood.val,ti.n_states);
    } else if(leaf_state[leaf_id]>ti.ambiguous_leaf_state) {
      die("invalid leaf state");

    } else { // regular state
      if(ti.tprob.is_branch_tprob_identity(branch_id)){
        site_prob[site_id]=prior_combined_lhood.val[leaf_state[leaf_id]];
      } else {
        const prob_t* const tprob_vec(branch_tprob+leaf_state[leaf_id]*ti.n_states);
        site_prob[site_id]=array_dot(prior_combined_lhood.val,tprob_vec,ti.n_states);
      }
    }

#ifdef LSC_DEBUG
    // check against conventional site lhood calculation:
    const prob_t check_prob(get_site_prob_single_cat(leaf_state,ti.tree,ti.tprob,ti.n_states));
    if(std::fabs(check_prob-site_cat_prob[site_id])>0.0001){
      log_os << "Error in site cached lhood calc: cached / norm: "
                << site_cat_prob[site_id] << " " << check_prob << "\n";
      abort();
    }
#endif
    // advance to the next unmasked site_id: must advance as far as
    // possible, even after break conditions are met:
    bool is_break(false);
    do {
      site_id++;
      if(site_id>=static_cast<int>(sdf.len)){
        is_break=true;
        break;
      }
      if(sdf.is_block_boundary[site_id]) is_break=true;
    } while(site_prob_mask[site_id]);

    if(is_break) break;
  }
}




static
void
eval_leaf(const tree_lhood_info& ti,
          const site_data_fastlup_core& sdf,
          const std::vector<tree_lhood_cache>& prior_branch_lhood,
          const int leaf_index,
          site_data_fastlup_core::index_type* const leaf_state,
          int& site_id,
          const std::vector<bool>& site_prob_mask,
          prob_t* const site_prob){

  const unsigned leaf_id(sdf.order[leaf_index]);
  const bi_tree_node* leaf_node(ti.tree.leaf_node(leaf_id));
  const bi_tree_node* parent_node(leaf_node->parent());
  const int branch_id(leaf_node->branch_id());
  const prob_t* const branch_tprob(ti.tprob.branch_prob(branch_id));

  if(parent_node==0) { die("Can't handle singleton tree"); }
  assert(leaf_index+1 < static_cast<int>(ti.tree.leaf_size()));

#ifdef LSC_DEBUG
  log_os << "Starting eval_leaf func; leaf_index, site_id " << leaf_index << " " << site_id << "\n";
  log_os << "leaf_node: " << leaf_node << "\n";
  log_os << "parent_node: " << parent_node << "\n";
#endif

  tree_walk_cache tw(ti.n_states);
  std::vector<tree_lhood_cache> local_prior_branch_lhood;

  while(site_id<static_cast<int>(sdf.len)){
    const site_data_fastlup_core::index_type* site_index(sdf.data_core[site_id].index);

    // find out if the state value of any previously evaluated leaf has changed at this site: if so
    // unwind the leaf-stack back down to that level:
    for(int i(leaf_index-1);i>=0;--i){
      const unsigned i_leaf_id(sdf.order[i]);
      if(leaf_state[i_leaf_id]!=site_index[i_leaf_id]) return;
    }

    leaf_state[leaf_id] = site_index[leaf_id];

#ifdef LSC_DEBUG
    log_os << "Site id loop; site_id: " << site_id << "\n";
    log_os << "leaf_index,leaf_state: " << leaf_index << " " << leaf_state[leaf_id] << "\n";
#endif

    // default branch_lhood for this org is a pointer into tprob
    tw.branch_lhood.node=leaf_node;

    if(leaf_state[leaf_id]==ti.ambiguous_leaf_state) { // ambiguous state case
      tw.branch_lhood.val=tw.branch_lhood_alloced;
      for(int i(0);i<ti.n_states;++i){
        tw.branch_lhood_alloced[i] = 1.;
      }

    } else if(leaf_state[leaf_id]>ti.ambiguous_leaf_state) {
      die("invalid leaf state");

    } else { // regular state case
      const prob_t* const leaf_branch_lhood(branch_tprob+leaf_state[leaf_id]*ti.n_states);
      tw.branch_lhood.val=leaf_branch_lhood;
    }

    // walk through tree to cache as much information as possible conditioned on this leaf state:
    local_prior_branch_lhood=prior_branch_lhood;
    lhood_site_cached_combine_branches(ti,parent_node,local_prior_branch_lhood,tw);
#ifdef LSC_DEBUG
    log_os << "tree walk returns branch/combined lhood node: " << tw.branch_lhood.node << " " << tw.combined_lhood.node << "\n";
    for(unsigned i(0);i<4;++i){ log_os << "branch_lhood: " << i << " " << tw.branch_lhood.val[i] << "\n"; }
#endif
    // with branch_lhood found, continue on to the state info for the next leaf:
    const unsigned next_leaf_index(leaf_index+1);
    if(next_leaf_index+1 == ti.tree.leaf_size()){
      eval_last_leaf(ti,sdf,tw.combined_lhood,next_leaf_index,leaf_state,site_id,site_prob_mask,site_prob);
    } else {
      local_prior_branch_lhood.push_back(tw.branch_lhood);
      eval_leaf(ti,sdf,local_prior_branch_lhood,next_leaf_index,leaf_state,site_id,site_prob_mask,site_prob);
    }
  }
}




//entry point function:
void
lhood_site_cached(const tree_lhood_info& ti,
                  const site_data_fastlup_core& sdf,
                  const std::vector<bool>& site_prob_mask,
                  prob_t* const site_prob){

  assert(ti.tree.leaf_size() == sdf.n_orgs);

  site_data_fastlup_core::index_type* leaf_state(new site_data_fastlup_core::index_type[ti.tree.leaf_size()]);

  static const int leaf_index(0);
  int site_id(0);

  // advance site_id to first unmasked position:
  while(site_id<static_cast<int>(sdf.len) && site_prob_mask[site_id]){ site_id++; }

  std::vector<tree_lhood_cache> prior_branch_lhood;

  eval_leaf(ti,sdf,prior_branch_lhood,leaf_index,leaf_state,site_id,site_prob_mask,site_prob);

  delete [] leaf_state;
}

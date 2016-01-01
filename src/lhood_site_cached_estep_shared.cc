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

// $Id: lhood_site_cached_estep_shared.cc 918 2007-10-12 23:36:13Z ctsa $
//

/// \file

#include "lhood_site_cached_estep_shared.h"
#include "lhood_site_cached_shared.h"
#include "util/math/matrix_util.h"
#include "util/general/die.h"

#ifdef LSC_DEBUG
#include "lhood_site.h"
#include "util/general/log.h"

#include <ostream>
#endif

void
lhood_site_cached_estep(const tree_lhood_info& ti,
                        const site_data_fastlup& sdf,
                        const std::vector<tree_lhood_cache>& prior_branch_lhood,
                        const tree_lhood_cache& prior_combined_lhood,
                        const int leaf_count,
                        site_data_fastlup_core::index_type* leaf_state,
                        int& site_id,
                        tree_walk_cache_estep_base& twe,
                        const bool is_reverse){

  // the current algorithm is only guaranteed to work for 3 species
  // trees, b/c it just travels back down the leaf stack, which does
  // not normally produce all required intermediate site probabilities
  // for the estep. will have to fix this if em minimization becomes
  // useful again
  //
  // elegantly notify user with a cryptic runtime assert:
  die("lhood_site_cached_estep: EM method disabled... see notes in code");

  const unsigned leaf_id(sdf.order[leaf_count]);
  const bi_tree_node* leaf_node(ti.tree.leaf_node(leaf_id));
  const bi_tree_node* parent_node(leaf_node->parent());
  const int branch_id(leaf_node->branch_id());
  const prob_t* const branch_tprob(ti.tprob.branch_prob(branch_id));

  if(parent_node==0) { die("Can't handle singleton tree"); }

#ifdef LSC_DEBUG
  log_os << "Entering estep eval_leaf func; leaf_count, site_id, reverse "
            << leaf_count << " " << site_id << " " << is_reverse << "\n";
  log_os << "leaf_node: " << leaf_node << "\n";
  log_os << "parent_node: " << parent_node << "\n";
#endif

  if( (leaf_count+1 == static_cast<int>(ti.tree.leaf_size())) && (! is_reverse) ){

    if(parent_node!=prior_combined_lhood.node){
      die("failed assertion in lhood_site_cached_estep");
    }

    leaf_state[leaf_id] = sdf.data_core[site_id].index[leaf_id];
#ifdef LSC_DEBUG
    log_os << "eval site_id: " << site_id << "\n";
    for(unsigned i(0); i<ti.tree.n_leaves(); ++i){
      log_os << "os: " << i << " " << ti.tree.leaf[i]->label << " " << leaf_state[i] << "\n";
    }
#endif

    if(leaf_state[leaf_id]==ti.ambiguous_leaf_state){ // ambiguous state
      die("can't handle ambiguous states using EM");
      twe.site_prob = array_sum(prior_combined_lhood.val,ti.n_states);
    } else if(leaf_state[leaf_id]>ti.ambiguous_leaf_state) {
      die("invalid leaf state");
    } else {  // regular state
      const prob_t* const tprob_vec(branch_tprob+leaf_state[leaf_id]*ti.n_states);
      twe.site_prob = array_dot(prior_combined_lhood.val,tprob_vec,ti.n_states);
    }

#ifdef LSC_DEBUG
    // check against conventional site lhood calculation:
    const prob_t check_prob(get_site_prob_single_cat(leaf_state,ti.tree,ti.tprob,ti.n_states));
    if(std::fabs(check_prob-twe.site_prob)>0.0001){
      log_os << "Error in site cached lhood calc: cached / norm: "
                << twe.site_prob << " " << check_prob << "\n";
      abort();
    }
#endif

    // now start counting back down leaves and repeat process in reverse order to fill in twe
    site_data_fastlup_core::index_type* leaf_state_rev(new site_data_fastlup_core::index_type[ti.tree.leaf_size()]);

    std::vector<tree_lhood_cache> local_prior_branch_lhood;
    lhood_site_cached_estep(ti,sdf,local_prior_branch_lhood,tree_lhood_cache(),leaf_count,leaf_state_rev,site_id,twe,true);

    delete [] leaf_state_rev;

    leaf_state[leaf_id] = sdf.data_core[site_id].index[leaf_id];

    // now that the double leaf walk is complete, calculate the statistics for this site
    // and store them ?somewhere in twe?
    twe.increment_site_stats(ti,sdf,site_id,leaf_state);

    site_id++;
    return;
  } else if( leaf_count==0 && is_reverse){
    //phew! have all info available for etrans analysis
    /// \todo it might be cleaner to find a way to call increment_site_stats here, rather than
    /// at the top of the reverse leaf stack -- what's currently in place gets the job done though,
    /// so no messing with it for now.
    ///
#ifdef LSC_DEBUG
    log_os << "made it! twe loaded.\n";
#endif
    return;
  }

  tree_walk_cache tw(ti.n_states);
  std::vector<tree_lhood_cache> local_prior_branch_lhood;

  while(site_id<static_cast<int>(sdf.len)){

    const site_data_fastlup_core::index_type* site_index(sdf.data_core[site_id].index);

    if(! is_reverse){
      // find out any lower org<leaf_count have changed at this site: if so
      // unwind the org-stack back down to that level:
      for(int i(leaf_count-1);i>=0;--i){
        const unsigned i_leaf_id(sdf.order[i]);
        if(leaf_state[i_leaf_id]!=site_index[i_leaf_id]) { return; }
      }
    }

    leaf_state[leaf_id] = site_index[leaf_id];

#ifdef LSC_DEBUG
    log_os << "Site id loop; site_id: " << site_id << "\n";
    log_os << "leaf_count,leaf_state: " << leaf_count << " " << leaf_state[leaf_id] << "\n";
#endif

    // default branch_lhood for this org is a pointer into tprob
    tw.branch_lhood.node=leaf_node;

    if(leaf_state[leaf_id]==ti.ambiguous_leaf_state){ // ambiguous state
      die("can't handle ambiguous states using EM");
      tw.branch_lhood.val=tw.branch_lhood_alloced;
      for(int i(0);i<ti.n_states;++i){
        tw.branch_lhood_alloced[i] = 1.;
      }
    } else if(leaf_state[leaf_id]>ti.ambiguous_leaf_state) {
      die("invalid leaf state");

    } else { // regular state
      const prob_t* const leaf_branch_lhood(branch_tprob+leaf_state[leaf_id]*ti.n_states);
      tw.branch_lhood.val=leaf_branch_lhood;
    }

    for(int i(0);i<ti.n_states;++i){
      twe.up_branch_prob[branch_id][i] = tw.branch_lhood.val[i];
    }

    // walk through tree to cache as much information as possible conditioned on this org state:
    local_prior_branch_lhood=prior_branch_lhood;
    lhood_site_cached_combine_branches(ti,parent_node,local_prior_branch_lhood,tw,&twe);
#ifdef LSC_DEBUG
    log_os << "tree walk returns branch/combined lhood nodes: " << tw.branch_lhood.node << " " << tw.combined_lhood.node << "\n";
    for(unsigned i(0);i<4;++i){ log_os << "branch lhood: " << i << " " << tw.branch_lhood.val[i] << "\n"; }
#endif
    // with branch_lhood found, continue on to the state info for the next org:
    unsigned next_leaf_count(leaf_count+1);
    if(is_reverse) next_leaf_count=leaf_count-1;

    local_prior_branch_lhood.push_back(tw.branch_lhood);
    lhood_site_cached_estep(ti,sdf,
                            local_prior_branch_lhood,
                            tw.combined_lhood,
                            next_leaf_count,leaf_state,site_id,twe,is_reverse);

    if(is_reverse) break;
  }
}

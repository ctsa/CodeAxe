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

// $Id: lhood_site_cached_estep.cc 1176 2008-03-26 00:17:56Z ctsa $

/// \file

#if 0
#include "lhood_site_cached_estep.h"
#include "lhood_site_cached_estep_shared.h"

// general tree em e-step -- like lhood_site_cached solution but after
// calculating site probability traverse back trough tree branches in
// reverse order to calculate this sites contribution to the expected
// substitutions per branch
//
// note that this is broken for anything but 3 species trees at the
// moment
//

static
void
local_increment_site_stats(const tree_lhood_info& ti,
                           const site_data_fastlup& sdf,
                           const tree_walk_cache_estep_base& twe,
                           const int site_id,
                           const site_data_fastlup_core::index_type* seq_state,
                           estep_dat_t& edat){

  const smlfloat pnorm(static_cast<smlfloat>(sdf.data[site_id].unassigned_cat_count)*(1./twe.site_prob));

  // do the root first
  for(int i(0);i<ti.n_states;++i){
    edat.root[i] +=
      pnorm*ti.tprob.root_prob()[i]*
      twe.up_branch_prob[ti.tree.root()->child1()->branch_id()][i]*
      twe.up_branch_prob[ti.tree.root()->child2()->branch_id()][i];
  }

  // now go through all branches:
  const unsigned n_branch(ti.tree.branch_size());
  for(unsigned i(0);i<n_branch;++i){
    const bi_tree_node* branch_node(ti.tree.branch_node(i));
    const int parent_branch_id(branch_node->parent()->branch_id());
    const int sister_branch_id(branch_node->sister()->branch_id());

    const prob_t* down_prob(ti.tprob.root_prob());
    if(parent_branch_id>=0) down_prob=twe.down_branch_prob[parent_branch_id];

    const prob_t* sister_prob(twe.up_branch_prob[sister_branch_id]);

    // deal with leaf branches separately:
    const int leaf_id(branch_node->leaf_id());
    if(leaf_id >= 0){

      const unsigned leaf_state(seq_state[leaf_id]);
       for(int a(0);a<ti.n_states;++a){
        prob_t val(pnorm*twe.up_branch_prob[i][a]*sister_prob[a]);
        edat.trans[i][leaf_state+a*ti.n_states] += val*down_prob[a];
      }

      // deal with internal branches:
    } else {
      for(int b(0);b<ti.n_states;++b){
        const prob_t val=pnorm*
          twe.up_branch_prob[branch_node->child1()->branch_id()][b]*
          twe.up_branch_prob[branch_node->child2()->branch_id()][b];

        const prob_t* tprob_vec(ti.tprob.branch_prob(i)+b*ti.n_states);

        for(int a(0);a<ti.n_states;++a){
          edat.trans[i][b+a*ti.n_states] +=
            val*down_prob[a]*
            tprob_vec[a]*
            sister_prob[a];
        }
      }
    }
  }
}



struct tree_walk_cache_estep : public tree_walk_cache_estep_base {

  typedef tree_walk_cache_estep_base base_t;

  tree_walk_cache_estep(const unsigned b,
                        const unsigned s,
                        estep_dat_t& init_edat)
    : base_t(b,s), edat(init_edat) {}

  virtual
  void
  increment_site_stats(const tree_lhood_info& ti,
                       const site_data_fastlup& sdf,
                       const int site_id,
                       const site_data_fastlup_core::index_type* seq_state){
    local_increment_site_stats(ti,sdf,*this,site_id,seq_state,edat);
  }

  estep_dat_t& edat;
};



//entry point function:
void
lhood_site_cached_estep(const tree_lhood_info& ti,
                        const site_data_fastlup& sdf,
                        estep_dat_t& edat){

  site_data_fastlup_core::index_type* leaf_state(new site_data_fastlup_core::index_type[ti.tree.leaf_size()]);

  int leaf_count(0);
  int site_id(0);

  tree_walk_cache_estep twe(ti.tree.branch_size(),ti.n_states,edat);

  std::vector<tree_lhood_cache> prior_branch_lhood;
  lhood_site_cached_estep(ti,sdf,prior_branch_lhood,tree_lhood_cache(),leaf_count,leaf_state,site_id,twe);

  delete [] leaf_state;
}
#endif

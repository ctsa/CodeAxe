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

// $Id: lhood_site_cached_estep_pprob.cc 918 2007-10-12 23:36:13Z ctsa $

/// \file

#include "lhood_site_cached_estep_pprob.h"
#include "lhood_site_cached_estep_shared.h"
#include "pprob_site.h"


/// elaborate estep_base with posterior prob values
///
struct tree_walk_cache_estep_pprob : public tree_walk_cache_estep_base {

  typedef tree_walk_cache_estep_base base_t;

  tree_walk_cache_estep_pprob(const unsigned n_branches,
                              const unsigned n_states,
                              pprob_cache& init_ppc)
    : base_t(n_branches,n_states), ppc(init_ppc) { }


  virtual
  void
  increment_site_stats(const tree_lhood_info& ti,
                       const site_data_fastlup&,
                       const int site_id,
                       const site_data_fastlup_core::index_type* seq_state){

    get_site_pprob_from_branch_prob(ti.tree,ti.n_states,ti.tprob.root_prob(),
                                    up_branch_prob.ptr(),down_branch_prob.ptr(),
                                    seq_state,ppc[site_id]);
  }

  pprob_cache& ppc;
};




//entry point function:
void
lhood_site_cached_estep_pprob(const tree_lhood_info& ti,
                              const site_data_fastlup& sdf,
                              pprob_cache& ppc){

  site_data_fastlup_core::index_type* seq_state(new site_data_fastlup_core::index_type[ti.tree.leaf_size()]);

  int org_id(0);
  int site_id(0);

  tree_walk_cache_estep_pprob twe(ti.tree.branch_size(),ti.n_states,ppc);

#ifdef DEBUG
  // check twe initialization w/ bogus fill-in:
  const unsigned vi(ti.tree.branch_size()*ti.n_states);
  std::fill(twe.up_branch_prob[0],twe.up_branch_prob[0]+vi,-1);
  std::fill(twe.down_branch_prob[0],twe.down_branch_prob[0]+vi,-1);
#endif

  std::vector<tree_lhood_cache> prior_branch_lhood;
  lhood_site_cached_estep(ti,sdf,prior_branch_lhood,tree_lhood_cache(),org_id,seq_state,site_id,twe);

  delete [] seq_state;
}

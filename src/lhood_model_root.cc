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

// $Id: lhood_model_root.cc 1031 2007-12-03 19:18:58Z ctsa $

/// \file

#include "cat_manager.h"
#include "condition_func.h"
#include "lhood_model_root.h"
#include "lhood_site.h"
#include "root_gtor.h"
#include "site_data_fastlup.h"
#include "subs_ml_model.h"
#include "util/math/array_util.h"
#include "util/math/matrix_util.h"



void
get_root_node_site_partial_prob_cat(const subs_ml_model& mdl,
                                    const site_data_fastlup& sdf,
                                    const unsigned cat_no,
                                    prob_t** root_spp,
                                    condition_func& cf){

  const unsigned n_states(mdl.state_size());
  const unsigned n_branches(mdl.tree().branch_size());
  const unsigned n_sites(sdf.len);

  simple_matrix<prob_t> dp_tmp(n_branches,n_states);

  tree_probs tprob;
  tprob.init_probs(mdl,cat_no);
  cf.tprob_init(tprob);
  for(unsigned i(0);i<n_sites;++i){
    get_site_root_partial_prob_single_cat_prep(sdf.data_core[i].index,
                                               mdl.tree(),
                                               tprob,
                                               n_states,
                                               root_spp[i],
                                               dp_tmp.ptr());
  }
}



void
get_root_node_site_partial_prob(const subs_ml_model& mdl,
                                const site_data_fastlup& sdf,
                                simple_matrix3d<prob_t>& root_spp,
                                std::auto_ptr<condition_func>* cf){

  const unsigned n_cats(mdl.get_cat_manager().cat_size());

  for(unsigned c(0);c<n_cats;++c){
    get_root_node_site_partial_prob_cat(mdl,sdf,c,root_spp[c],*cf[c]);
  }
}



void
get_site_prob_from_root_spp(const root_gtor& rg,
                            const prob_t * const * root_spp,
                            const unsigned n_sites,
                            const unsigned cat_no,
                            prob_t* site_prob){

  const unsigned n_states(rg.state_size());

  const prob_t* root_p(rg.cat_state_pdistro(cat_no));

  std::fill(site_prob,site_prob+n_sites,0.);
  transposed_matrix_vector_mult_sum(root_spp[0],root_p,site_prob,n_sites,n_states);
}

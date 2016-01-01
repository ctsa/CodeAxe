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

// $Id: lhood_model_prep_submodel.cc 1222 2008-05-22 23:10:06Z ctsa $

/// \file

#include "bi_tree.h"
#include "cat_manager.h"
#include "condition_func.h"
#include "lhood_model_prep_submodel.h"
#include "rate_gtor.h"
#include "site_data_fastlup.h"
#include "site_data_fastlup_core_util.h"
#include "site_data_fastlup_util.h"
#include "subs_ml_model.h"
#include "tree_probs.h"
#include "util/general/die.h"



void
submodel_data::
init(const subs_ml_model& mdl,
     const int init_n_states,
     const site_data_fastlup& full_sdf,
     const unsigned* index_translation){

  if(is_init){ die("submodel_data bad init"); }

  const bi_tree& tree(mdl.tree());

  n_states = init_n_states;
  full_sdf_index_map.init(full_sdf.len);

  {
    site_data_fastlup tmp_sdf;
    site_data_fastlup_state_reduction(tree,n_states,full_sdf,index_translation,tmp_sdf);
    sdf = tmp_sdf;

    // create cat masks from temporary ads_masks:
    const cat_manager& cm(mdl.get_cat_manager());
    const unsigned n_cats(cm.cat_size());
    const unsigned n_adsets(cm.assigned_data_set_size());

    simple_matrix<bool> is_adset_using_cat(n_adsets,n_cats);
    for(unsigned a(0);a<n_adsets;++a){
      cm.is_adset_using_cat(is_adset_using_cat[a],a);
    }

    cat_site_mask.resize(n_cats);

    for(unsigned c(0);c<n_cats;++c){
      cat_site_mask[c].resize(tmp_sdf.len,true);
      for(unsigned a(0);a<n_adsets;++a){
        if(is_adset_using_cat[a][c]) {
          std::transform(cat_site_mask[c].begin(),cat_site_mask[c].end(),
                         tmp_sdf.ads_site_mask(a).begin(),
                         cat_site_mask[c].begin(),std::logical_and<bool>());
        }
      }
    }
  }

  map_translated_index_superset_data(full_sdf,sdf,index_translation,full_sdf_index_map.ptr());

  is_init=true;
}



submodel_manager::
submodel_manager(const subs_ml_model& mdl,
                 const site_data_fastlup& full_sdf)
  : n_submodels(static_cast<int>(mdl.submodel_size())),
    site_prob(n_submodels),
    _submodel(n_submodels),_tprob(n_submodels),
    _cf(n_submodels){

  const int n_states(static_cast<int>(mdl.state_size()));
  simple_array<unsigned> index_translation(n_states+1);

  for(int s(0);s<n_submodels;++s){
    mdl.get_rate_gtor().submodel_state_reduction_map(s,index_translation.ptr());

    const int submodel_n_states(static_cast<int>(mdl.submodel_state_size(s)));
    submodel(s).init(mdl,submodel_n_states,full_sdf,index_translation.ptr());
    site_prob[s].init(submodel(s).sdf.len);
    _cf[s].reset(mdl.get_rate_gtor().submodel_condition_func_factory(s));
    cf(s).data_init(mdl.tree(),submodel(s).sdf);
  }
}

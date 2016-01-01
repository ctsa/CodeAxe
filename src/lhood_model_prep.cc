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

// $Id: lhood_model_prep.cc 1222 2008-05-22 23:10:06Z ctsa $

/// \file

#include "cat_manager.h"
#include "condition_func.h"
#include "lhood_model_prep.h"
#include "lhood_model_prep_group.h"
#include "lhood_model_prep_submodel.h"
#include "rate_gtor.h"
#include "site_data_fastlup.h"
#include "subs_ml_model.h"
#include "tree_probs.h"



lhood_model_prep_cat_site_prob::
lhood_model_prep_cat_site_prob(const subs_ml_model& mdl,
                               const site_data_fastlup& sdf)
  : _tprob(new tree_probs), site_prob(0), tprob_ws(3) {

  const unsigned full_count_len(sdf.len);
  const bool is_use_submodels(mdl.submodel_size()>1);

  if(is_use_submodels){
    _sm.reset(new submodel_manager(mdl,sdf));
  } else {
    _cf.reset(mdl.get_rate_gtor().condition_func_factory());
    cf().data_init(mdl.tree(),sdf);
  }

  site_prob = new prob_t[full_count_len];

  // prepare category specific masks by bit-anding together all
  // ads masks for each category:
  //
  const cat_manager& cm(mdl.get_cat_manager());
  const unsigned n_cats(cm.cat_size());
  const unsigned n_adsets(cm.assigned_data_set_size());

  // setup is_a_u_c:
  is_adset_using_cat.init(n_adsets,n_cats);

  for(unsigned a(0);a<n_adsets;++a){
    cm.is_adset_using_cat(is_adset_using_cat[a],a);
  }

  // setup cat_site_mask
  cat_site_mask.resize(n_cats);

  for(unsigned c(0);c<n_cats;++c){
    cat_site_mask[c].resize(full_count_len,true);
    for(unsigned a(0);a<n_adsets;++a){
      if(is_adset_using_cat[a][c]) {
        std::transform(cat_site_mask[c].begin(),cat_site_mask[c].end(),
                       sdf.ads_site_mask(a).begin(),
                       cat_site_mask[c].begin(),std::logical_and<bool>());
      }
    }
  }
}



lhood_model_prep_cat_site_prob::
~lhood_model_prep_cat_site_prob(){
  delete [] site_prob;
}



lhood_model_prep::
lhood_model_prep(const subs_ml_model& mdl,
                 const site_data_fastlup& sdf)
  : csp_prep(mdl,sdf) {

  typedef site_data_fastlup_cell::group_assigned_data_set_vector::const_iterator gavi;

  const cat_manager& cm(mdl.get_cat_manager());
  const unsigned n_cats(cm.cat_size());
  const unsigned n_group_cats(cm.group_cat_size());
  const unsigned n_adsets(cm.assigned_data_set_size());
  const unsigned n_groups(sdf.n_groups);

  _gd.reset(new group_data(n_groups));

  site_cat_mix_prob.init(n_adsets,sdf.len);

  cat_pdistro.resize(n_cats);
  group_cat_pdistro.resize(n_group_cats);

  adset_cat_pdistro.init(n_adsets,n_cats);
  adset_group_cat_pdistro.init(n_adsets,n_group_cats);

  // with multiple adsets, group priors can become dependend on
  // group_id whenever a group contains sites assigned to different
  // adsets
  //
  _prob_adset_on_group.init(n_groups,n_adsets,0.);

  for(unsigned i(0);i<sdf.len;++i){
    gavi g(sdf.data[i].group_assigned_data_set_count.begin());
    const gavi g_end(sdf.data[i].group_assigned_data_set_count.end());
    for(;g!=g_end;++g){
      const site_data_fastlup_code::group_id_type group_id(g->first.group_id);
      const site_data_fastlup_code::assigned_data_set_id_type adset_id(g->first.assigned_data_set_id);
      const unsigned count(g->second);

      _prob_adset_on_group[group_id][adset_id] += count;
    }
  }

  for(unsigned g(0);g<n_groups;++g){
    pdistro_norm(_prob_adset_on_group[g],_prob_adset_on_group[g]+n_adsets);
  }

  ads_group_cat_prior_norm.resize(n_adsets);
  group_id_group_cat_pdistro.init(n_groups,n_group_cats);
}



// leave empty destructor here so that auto_ptr's can close out correctly:
lhood_model_prep::~lhood_model_prep(){}

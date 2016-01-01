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

// $Id: site_prob_maker.cc 1175 2008-03-26 00:04:19Z ctsa $

/// \file

#include "condition_func.h"
#include "rate_gtor.h"
#include "lhood_model_prep.h"
#include "lhood_model_prep_submodel.h"
#include "lhood_model_root.h"
#include "lhood_site_cached.h"
#include "site_data_fastlup.h"
#include "site_prob_maker.h"
#include "subs_ml_model.h"
#include "subs_ml_model_min_options.h"
#include "tree_probs.h"



/// \brief convenience wrapper for lhood_site_cached call
///
static
void
get_cached_site_prob_array(const bi_tree& tree,
                           const tree_probs& tprob,
                           const site_data_fastlup_core& sdf,
                           const int n_states,
                           const condition_func& cf,
                           const std::vector<bool>& site_prob_mask,
                           prob_t* const site_prob){

  const tree_lhood_info ti(tree,tprob,n_states);

  lhood_site_cached(ti,sdf,site_prob_mask,site_prob);

  cf.condition_site_prob(tprob,site_prob_mask,site_prob);
}



void
site_prob_maker_full::
get_cat_site_prob(const subs_ml_model& mdl,
                  const unsigned cat_no,
                  lhood_model_prep_cat_site_prob& nlpc) const {

  const bool is_use_submodels(mdl.submodel_size()>1);
  const rate_gtor& rg(mdl.get_rate_gtor());

  std::fill(nlpc.site_prob,nlpc.site_prob+_full_sdf.len,1.);

  if(! is_use_submodels) {
    nlpc.tprob().init_probs(mdl,cat_no,is_use_submodels,0,nlpc.tprob_ws.ptr());
    get_cached_site_prob_array(mdl.tree(),nlpc.tprob(),_full_sdf,
                               mdl.state_size(),nlpc.cf(),nlpc.cat_site_mask[cat_no],
                               nlpc.site_prob);
  } else {

    for(int s(0);s<nlpc.sm().n_submodels;++s){
      std::fill(nlpc.sm().site_prob[s].ptr(),nlpc.sm().site_prob[s].ptr()+nlpc.sm().submodel(s).sdf.len,0.);
      nlpc.sm().tprob(s).init_probs(mdl,cat_no,is_use_submodels,s,nlpc.tprob_ws.ptr());

      const submodel_data& sd(nlpc.sm().submodel(s));
      get_cached_site_prob_array(mdl.tree(),nlpc.sm().tprob(s),sd.sdf,
                                 sd.n_states,nlpc.sm().cf(s),sd.cat_site_mask[cat_no],
                                 nlpc.sm().site_prob[s].ptr());

      for(unsigned i(0);i<_full_sdf.len;++i){
        const int s_index(sd.full_sdf_index_map[i]);
        nlpc.site_prob[i] *= nlpc.sm().site_prob[s][s_index];
      }
    }
    rg.submodel_adjust_p(nlpc.site_prob,_full_sdf.len);
  }
}



void
site_prob_maker_root::
get_cat_site_prob(const subs_ml_model& mdl,
                  const unsigned cat_no,
                  lhood_model_prep_cat_site_prob& nlpc) const {

  const bool is_use_submodels(mdl.submodel_size()>1);
  if(is_use_submodels) pass_away("submodel approx not supported by root-cycler");

  const unsigned n_sites(_root_spp.dim2());
  get_site_prob_from_root_spp(mdl.get_root_gtor(),_root_spp[cat_no],n_sites,cat_no,nlpc.site_prob);
  _cf[cat_no]->condition_site_prob_prep(nlpc.cat_site_mask[cat_no],nlpc.site_prob);
}



void
site_prob_maker_rootcat::
get_cat_site_prob(const subs_ml_model& mdl,
                  const unsigned cat_no,
                  lhood_model_prep_cat_site_prob& nlpc) const {

  const bool is_use_submodels(mdl.submodel_size()>1);
  if(is_use_submodels) pass_away("submodel approx not supported by root-cycler");

  const unsigned n_sites(_root_spp.dim1());

  if(cat_no != _cat_no){
    const unsigned ocat_no( (cat_no<_cat_no) ? cat_no : cat_no-1);
    std::copy(_osp[ocat_no],_osp[ocat_no]+n_sites,nlpc.site_prob);
  } else {
    get_site_prob_from_root_spp(mdl.get_root_gtor(),_root_spp.ptr(),n_sites,cat_no,nlpc.site_prob);
    _cf.condition_site_prob_prep(nlpc.cat_site_mask[cat_no],nlpc.site_prob);
  }
}

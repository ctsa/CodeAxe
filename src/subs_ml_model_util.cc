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

// $Id: subs_ml_model_util.cc 1192 2008-03-28 23:41:48Z ctsa $

/// \file

#include "bg_gtor.h"
#include "bi_tree.h"
#include "cat_info.h"
#include "cat_manager.h"
#include "data_class_assignment_util.h"
#include "rate_gtor_nuc_base.h"
#include "rate_gtor_nscodon_base.h"
#include "rate_gtor_util.h"
#include "root_gtor.h"
#include "simple_util.h"
#include "site_data.h"
#include "site_data_stats.h"
#include "subs_ml_model.h"
#include "subs_ml_ptol.h"


void
subs_ml_model::
init_trained_data_dependencies(const site_data& sd){

  const cat_manager& cm(get_cat_manager());
  const unsigned n_states(state_size());
  const unsigned n_adsets(cm.assigned_data_set_size());

  if(n_adsets==1){
    simple_array<prob_t> distro(n_states);
    site_data_2_state_pdf(sd,n_states,distro.ptr());

    get_root_gtor_nonconst().set_pdistro(distro.ptr());
    get_bg_gtor_nonconst().set_pdistro(distro.ptr());

  } else {

    die("not updated for bg_gtor -- and root gtor update doesn't make sense with assigned data sets");

#if 0
    const unsigned n_cats(cm.cat_size());
    const unsigned n_root_cats(cm.typed_cat_size(CAT_PARAM_TYPE::ROOT));
    if(n_adsets != n_cats) die("assigned cat/ model cat mismatch");
    if(n_adsets != n_root_cats)
      die("model cat -- root cat mismatch.. write more code to handle this!!");

    const unsigned n_data_classes(sd.data_class_size());

    data_class_id_assignment_map dc_ads_map(cm,sd);

    //
    // set the root distros from the corresponding assigned_data_set distros
    //
    std::vector<unsigned> data_class_state_count;
    site_data_2_data_class_state_count(sd,n_data_classes,n_states,data_class_state_count);

    simple_matrix<prob_t> assigned_data_set_state_distro(n_adsets,n_states,0.);

    for(unsigned i(0);i<n_data_classes;++i){
      if(! dc_ads_map.is_data_class_assigned(i)) continue;
      const unsigned adset_id(dc_ads_map.get_data_class_assigned_data_set_id(i));

      for(unsigned j(0);j<n_states;++j){
        assigned_data_set_state_distro[adset_id][j] += data_class_state_count[i*n_states+j];
      }
    }

    // finish this section to map root cats to adsets somehow:
    const unsigned n_root_cats(XXX);
    for(unsigned i(0);i<n_adsets;++i){
      pdistro_norm(model_cat_state_distro[i],model_cat_state_distro[i]+model_cat_state_distro.dim2());
      get_root_gtor_nonconst().set_distro_cat_state_pdistro(i,model_cat_state_distro[i]);
    }
#endif
  }
}



void
subs_ml_model::
init_untrained_data_dependencies(const site_data& sd){

  const cat_manager& cm(get_cat_manager());
  const unsigned n_states(state_size());
  const unsigned n_data_classes(sd.data_class_size());
  const unsigned n_adsets(cm.assigned_data_set_size());
#if 0
  std::vector<unsigned> data_class_state_count;
  site_data_2_data_class_state_count(sd,n_data_classes,n_states,data_class_state_count);
#endif

  data_class_id_assignment_map dc_ads_map(cm,sd);

  const unsigned n_seqs(sd.n_taxa());

  simple_matrix3d<unsigned> adset_seq_state_counts(n_adsets,n_seqs,n_states,0);

  for(unsigned i(0);i<n_data_classes;++i){
    if(! dc_ads_map.is_data_class_assigned(i)) continue;
    const unsigned adset_id(dc_ads_map.get_data_class_assigned_data_set_id(i));
    site_data_2_taxa_state_count(sd,n_states,adset_seq_state_counts[adset_id][0],true,true,i);
  }

  // map from sd seqs to tree leaves:
  std::vector<unsigned> sd_org_tree_order(n_seqs);
  for(unsigned j(0);j<n_seqs;++j){
    sd_org_tree_order[j] = sd.taxid.getid(tree().leaf_node(j)->label());
  }

  // remap sequence order
  simple_matrix<unsigned> ptmp(n_seqs,n_states);
  for(unsigned i(0);i<n_adsets;++i){
    unsigned** si(adset_seq_state_counts[i]);
    for(unsigned j(0);j<n_seqs;++j){
      for(unsigned k(0);k<n_states;++k){
        ptmp[j][k] = si[j][k];
      }
    }
    for(unsigned j(0);j<n_seqs;++j){
      for(unsigned k(0);k<n_states;++k){
        si[j][k] = ptmp[sd_org_tree_order[j]][k];
      }
    }
  }

  // remap assigned cats distros to obs cats:
  const unsigned n_obs_cats(cm.typed_cat_size(CAT_PARAM_TYPE::OBS));
  if(n_adsets==1){
    simple_array<prob_t> obs_cp(n_obs_cats);
    cm.typed_cat_pdistro(obs_cp.ptr(),CAT_PARAM_TYPE::OBS);

    simple_matrix3d<smlfloat> ssc(n_obs_cats,n_seqs,n_states);
    unsigned** s(adset_seq_state_counts[0]);
    for(unsigned i(0);i<n_obs_cats;++i){
      const prob_t p(obs_cp[i]);
      for(unsigned j(0);j<n_seqs;++j){
        for(unsigned k(0);k<n_states;++k){
          ssc[i][j][k] = p*static_cast<smlfloat>(s[j][k]);
        }
      }
    }

    set_obs_cat_seq_state_counts(ssc.ptr());
  } else {
    if(n_obs_cats != n_adsets){
      die("can't handle general obs cat<->assigned data set combinations");
    }
    simple_matrix3d<smlfloat> ssc(n_obs_cats,n_seqs,n_states);
    for(unsigned i(0);i<n_obs_cats;++i){
      unsigned** si(adset_seq_state_counts[i]);
      for(unsigned j(0);j<n_seqs;++j){
        for(unsigned k(0);k<n_states;++k){
          ssc[i][j][k] = static_cast<smlfloat>(si[j][k]);
        }
      }
    }
    set_obs_cat_seq_state_counts(ssc.ptr());
  }
}



void
subs_ml_model::
norm(){

  const cat_manager& cm(get_cat_manager());

  { // test no-go conditions:
    const unsigned n_cats(cm.cat_size());
    for(unsigned c(0);c<n_cats;++c){
      const unsigned n_branch_cat_sets(cm.branch_cat_set_size(c));
      if(n_branch_cat_sets>1){
        pass_away("model equilibration not possible with branch category models.");
      }
    }

    const unsigned n_adsets(cm.assigned_data_set_size());
    if(n_adsets > 1){
      /// \todo fix norm_model to work with assigned cats
      pass_away("model equilibration has not yet been updated for multiple assigned data sets");
    }
  }

#if 0
  simple_matrix<smlfloat> assigned_cat_distro(n_assigned_cats,n_states);
  for(unsigned i(0);i<n_assigned_cats;++i){
    array_zero(assigned_cat_distro.val[i],n_states);
  }
#endif

  static const unsigned assumed_branch_cat_set(0);

  const unsigned n_states(state_size());
  const rate_gtor& rg(get_rate_gtor());
  const root_gtor& ro(get_root_gtor());

  const unsigned n_root_distro_cats(ro.distro_cat_size());
  simple_matrix<prob_t> root_cat_state_pdistro(n_root_distro_cats,n_states,0.);

  const unsigned n_obs_cats(cm.typed_cat_size(CAT_PARAM_TYPE::OBS));
  simple_matrix<prob_t> obs_cat_state_pdistro(n_obs_cats,n_states,0.);

  const unsigned n_bg_distro_cats(get_bg_gtor().distro_cat_size());

  if(n_bg_distro_cats != n_obs_cats){
    die("parameterized bg_gtor not currently handled in norm()");
  }

  {
    simple_array<smlfloat> cat_state_pdistro(n_states);
    const unsigned n_cats(cm.cat_size());
    simple_array<prob_t> cp(n_cats);
    cm.cat_pdistro(cp.ptr());

    for(unsigned c(0);c<n_cats;++c){
      get_stationary_pdistro(cat_state_pdistro.ptr(),rg,
                             rates_func_options(rates_func_options_base(c),assumed_branch_cat_set));
      pdistro_check(cat_state_pdistro.ptr(),n_states,SUBS_ML_PTOL);

      const unsigned root_distro_cat(ro.get_distro_cat_no_from_cat_no(c));
      const unsigned obs_cat(cm.typed_cat_no_from_cat_no(c,CAT_PARAM_TYPE::OBS));

      const prob_t cat_prob(cp[c]);
      for(unsigned j(0);j<n_states;++j){
        const prob_t val(cat_state_pdistro[j]*cat_prob);
        obs_cat_state_pdistro[obs_cat][j] += val;
        root_cat_state_pdistro[root_distro_cat][j] += val;
      }
    }
  }

  root_gtor& rox(get_root_gtor_nonconst());
  for(unsigned i(0);i<n_root_distro_cats;++i){
    pdistro_norm(root_cat_state_pdistro[i],
                 root_cat_state_pdistro[i]+root_cat_state_pdistro.dim2());
    rox.set_distro_cat_state_pdistro(i,root_cat_state_pdistro[i]);
  }

  const unsigned n_seqs(tree().leaf_size());
  simple_matrix<prob_t*> oc(n_obs_cats,n_seqs);

  for(unsigned i(0);i<n_obs_cats;++i){
    pdistro_norm(obs_cat_state_pdistro[i],
                 obs_cat_state_pdistro[i]+obs_cat_state_pdistro.dim2());

    for(unsigned j(0);j<n_seqs;++j) oc[i][j]=obs_cat_state_pdistro[i];
  }
  set_obs_cat_seq_state_distro(oc.ptr());
}

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

// $Id: tree_probs.cc 1160 2008-03-20 19:18:02Z ctsa $

/// \file

#include "bi_tree.h"
#include "branch_time_util.h"
#include "root_gtor.h"
#include "subs_ml_model.h"
#include "subs_ml_ptol.h"
#include "subs_prob_util.h"
#include "subs_prob_util_core.h"
#include "time_gtor.h"
#include "tree_probs.h"
#include "util/math/matrix_util.h"
#include "util/math/prob_util_io.h"

#ifdef BRANCH_SPECIFIC_OBS
#include "rate_edge_dependencies.h"
#include "subs_ml_model_root_from_obs.h"
#endif


void
tree_probs::
alloc(const unsigned init_n_states,
      const unsigned init_n_branches){

  n_states=init_n_states;
  n_branches=init_n_branches;

  _root_prob.init(n_states);
  _branch_prob.init(n_branches,n_states*n_states);

  _is_branch_tprob_identity.init(n_branches,false);
}



void
tree_probs::
update_root(const prob_t* root){
  static const bool is_use_submodel(false);
  if(is_use_submodel){
    die("can't handle submodel root update");
  } else {
    std::copy(root,root+n_states,root_prob());
  }
}



void
tree_probs::
init_probs(const subs_ml_model& mdl,
           const unsigned cat,
           const bool is_use_submodel,
           const unsigned submodel_no,
           workspace<char> ws[3]){

  _is_site_prob_mode=false;

  const unsigned& init_n_states(is_use_submodel ? mdl.submodel_state_size(submodel_no):mdl.state_size());
  const unsigned& init_n_branches(mdl.tree().branch_size());

  if(n_states != init_n_states ||
     n_branches != init_n_branches){
    alloc(init_n_states,init_n_branches);
  }

  const prob_t* full_root(mdl.get_root_gtor().cat_state_pdistro(cat));
  if(is_use_submodel){
    mdl.get_rate_gtor().submodel_pdistro_reduction(full_root,submodel_no,root_prob());
  } else {
    std::copy(full_root,full_root+n_states,root_prob());
  }

#ifndef NDEBUG
  pdistro_check(root_prob(),n_states,SUBS_ML_PTOL);
#endif

  const rates_func_options_base ropt(cat,false,false,is_use_submodel,submodel_no);

#ifdef BRANCH_SPECIFIC_OBS
  die("branch specific obs not yet updated to account for branch categories");
  const prob_t * const * node_state_prob(mdl.get_root_gtor().cat_node_state_prob(cat));

  site_model_edge_pdistros smp;
  if(node_state_prob){
    smp.ropt_hookup(ropt);
  }

  const SITE_MODEL::index_t sm(mdl.get_rate_gtor().site_model());
  const time_gtor& tgm(mdl.get_time_gtor());

  if(node_state_prob){
    const prob_t * const * zzz(branch_state_prob);
    simple_array<prob_t> state_prob_tmp(n_states);
    simple_array<prob_t> tprob_tmp(n_states*n_states);
    simple_array<prob_t> tprob_tmp2(n_states*n_states);

    for(unsigned b(0);b<n_branches;++b){
      static const unsigned n_sections(3);
      const smlfloat section_frac(1./static_cast<smlfloat>(n_sections));
      const smlfloat section_time(tgm.branch_time(b)*section_frac);

      for(unsigned j(0);j<n_sections;++j){
        const unsigned node_id(b+1);
        const unsigned parent_node_id(mdl.tree().get_parent_index(node_id));
        const smlfloat node_frac((section_frac*0.5)*(1+j*2));

        for(unsigned i(0);i<n_states;++i){
          state_prob_tmp.val[i] = zzz[node_id][i]*(node_frac)+
            zzz[parent_node_id][i]*(1.-node_frac);
        }
        set_site_model_edge_pdistros(sm,state_prob_tmp.val,smp);
        get_single_branch_subs_prob(tprob_tmp.val,section_time,mdl,ropt,b,ws);

        if(j==0){
          matrix_copy(tprob_tmp.val,branch_prob(b),n_states);
        } else {
          matrix_mult(tprob_tmp.val,tprob_tmp2.val,branch_prob(b),n_states);
        }
        if((j+1)!=n_sections){
          matrix_copy(branch_prob(b),tprob_tmp2.val,n_states);
        }
      }
    }

  } else {
    for(unsigned b(0);b<n_branches;++b){
      get_single_branch_subs_prob(branch_prob(b),tgm.branch_time(b),mdl,ropt,b,ws);
    }
  }
#else

  get_all_branch_subs_prob(_branch_prob.ptr(),mdl,ropt,ws);

#endif

  for(unsigned b(0);b<n_branches;++b){
    _is_branch_tprob_identity[b]=is_zero_branch_time(mdl.get_time_gtor().branch_time(b,cat));
  }
  for(unsigned b(0);b<n_branches;++b){
    matrix_transpose_inplace(branch_prob(b),n_states);
  }
}



void
tree_probs::
self_check() const {

  pdistro_check(root_prob(),n_states,SUBS_ML_PTOL);
  for(unsigned b(0);b<n_branches;++b){
    for(unsigned i(0);i<n_states;++i){
      pdistro_check(branch_prob(b)+i,n_states,SUBS_ML_PTOL,n_states);
    }
  }
}

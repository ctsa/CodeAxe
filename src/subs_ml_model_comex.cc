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

// $Id: subs_ml_model_comex.cc 1153 2008-03-18 00:24:07Z ctsa $

/// \file

#include "bg_gtor.h"
#include "subs_ml_model.h"
#include "subs_ml_model_comex.h"
#include "subs_ml_model_root_from_obs.h"



const bi_tree&
cat_manager_sml_share::
tree() const { return _smlm.tree(); }



const cat_manager&
tr_gtor_sml_share::
get_cat_manager() const { return _smlm.get_cat_manager(); }



const bi_tree&
time_gtor_sml_share::
tree() const { return _smlm.tree(); }

bi_tree&
time_gtor_sml_share::
tree_nonconst() { return _smlm_nonconst.tree_nonconst(); }




RATE_GTOR_MODEL::index_t
r_gtor_sml_share::
rate_gtor_model() const { return _smlm.rate_gtor_model(); }

unsigned
r_gtor_sml_share::
state_size() const { return _smlm.state_size(); }



const notifier&
rate_gtor_sml_share::
get_bg_gtor_notifier() const {
  return _smlm.bg_gtor_notifier();
}

unsigned
rate_gtor_sml_share::
bg_cat_size() const{
  return _smlm.get_bg_gtor().distro_cat_size();
}

unsigned
rate_gtor_sml_share::
bg_cat_no_from_cat_no(const unsigned cat_no) const {
  return _smlm.get_bg_gtor().get_distro_cat_no_from_cat_no(cat_no);
}

const prob_t*
rate_gtor_sml_share::
bg_state_pdistro_bg_cat(const unsigned bg_cat_no) const {
  return _smlm.get_bg_gtor().distro_cat_state_pdistro(bg_cat_no);
}



unsigned
prob_gtor_obs_sml_share::
state_size() const { return _smlm.state_size(); }


void
prob_gtor_obs_sml_share::
get_obs_state_avg_root(prob_t* root_down_prob,
                       const unsigned cat_no) const {

  const prob_t* obs(_smlm.obs_state_distro_cat(cat_no));

  /// \todo: this is a gigantic waste: b/c we could just take the const*,
  /// how can the prob_gtor that uses this be speciallized?
  std::copy(obs,obs+state_size(),root_down_prob);
}

void
prob_gtor_obs_sml_share::
get_obs_state_time_avg_root(prob_t* root_down_prob,
                            const unsigned cat_no) const {
  subs_ml_model_obs_state_time_avg_root(_smlm,root_down_prob,cat_no);
}

bool
prob_gtor_obs_sml_share::
get_root_node_prob(prob_t* root_node_prob,
                   const unsigned cat_no) const {
  return subs_ml_model_root_node_prob(_smlm,root_node_prob,cat_no);
}

#ifdef BRANCH_SPECIFIC_OBS
unsigned
prob_gtor_obs_sml_share::
branch_size() const {
  return _smlm.tree().n_branches();
}

bool
prob_gtor_obs_sml_share::
get_root_node_prob_nspround(prob_t* root_node_prob,
                            const unsigned cat_no,
                            prob_t** node_state_prob) const {

  return subs_ml_model_root_node_prob_nspround(_smlm,root_node_prob,
                                               cat_no,node_state_prob);
}
#endif

const notifier&
prob_gtor_obs_sml_share::
get_obs_info_notifier() const { return _smlm.obs_info_notifier(); }



const prob_t*
root_gtor_sml_share::
obs_state_distro_cat(const unsigned cat_no) const { return _smlm.obs_state_distro_cat(cat_no); }



const notifier&
bg_gtor_sml_share::
get_obs_info_notifier() const { return _smlm.obs_info_notifier(); }

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

// $Id: subs_ml_model_comex.h 1153 2008-03-18 00:24:07Z ctsa $

/// \file

#ifndef __SUBS_ML_MODEL_COMEX_H
#define __SUBS_ML_MODEL_COMEX_H

#include "rate_gtor_options.h"
#include "subs_ml_types.h"

struct bi_tree;
struct cat_manager;
struct notifier;
struct subs_ml_model;
struct time_gtor;


/// \brief subs_ml_model component exchange
///
/// present a restricted subset of subs_ml_model functions to model
/// components, such that component communication must move through
/// this structure
///
struct subs_ml_model_comex {
  typedef subs_ml_model_comex self_t;

  subs_ml_model_comex(const subs_ml_model& smlm) : _smlm(smlm) {}

private:
  self_t& operator=(const self_t&);

protected:
  const subs_ml_model& _smlm;
};


struct cat_manager_sml_share : public subs_ml_model_comex {
  typedef subs_ml_model_comex base_t;

  cat_manager_sml_share(const subs_ml_model& smlm) : base_t(smlm) {}

  const bi_tree& tree() const;
};


struct tr_gtor_sml_share : public subs_ml_model_comex {
  typedef subs_ml_model_comex base_t;

  tr_gtor_sml_share(const subs_ml_model& smlm) : base_t(smlm) {}

  const cat_manager& get_cat_manager() const;
};


struct time_gtor_sml_share : public tr_gtor_sml_share {
  typedef tr_gtor_sml_share base_t;

  time_gtor_sml_share(subs_ml_model& smlm) : base_t(smlm), _smlm_nonconst(smlm) {}

  bi_tree& tree_nonconst();
  const bi_tree& tree() const;

private:
  subs_ml_model& _smlm_nonconst;
};


struct r_gtor_sml_share : public tr_gtor_sml_share {
  typedef tr_gtor_sml_share base_t;

  r_gtor_sml_share(const subs_ml_model& smlm) : base_t(smlm) {}

  RATE_GTOR_MODEL::index_t rate_gtor_model() const;
  unsigned state_size() const;
};

struct rate_gtor_sml_share : public r_gtor_sml_share {
  typedef r_gtor_sml_share base_t;

  rate_gtor_sml_share(const subs_ml_model& smlm) : base_t(smlm) {}

  const notifier& get_bg_gtor_notifier() const;

  unsigned bg_cat_size() const;

  unsigned bg_cat_no_from_cat_no(const unsigned cat_no) const;

  const prob_t* bg_state_pdistro_bg_cat(const unsigned bg_cat_no) const;
};



struct prob_gtor_obs_sml_share : public subs_ml_model_comex {
  typedef subs_ml_model_comex base_t;

  prob_gtor_obs_sml_share(const subs_ml_model& smlm) : base_t(smlm) {}

  const notifier& get_obs_info_notifier() const;

  unsigned state_size() const;

  void get_obs_state_avg_root(prob_t* root_down_prob,
                              const unsigned cat_no) const;

  void get_obs_state_time_avg_root(prob_t* root_down_prob,
                                   const unsigned cat_no) const;

  bool get_root_node_prob(prob_t* root_node_prob,
                          const unsigned cat_no) const;

#ifdef BRANCH_SPECIFIC_OBS
  unsigned branch_size() const;

  bool get_root_node_prob_nspround(prob_t* root_node_prob,
                                   const unsigned cat_no,
                                   prob_t** node_state_prob) const;
#endif
};



struct multi_prob_gtor_sml_share : public r_gtor_sml_share {
  typedef r_gtor_sml_share base_t;

  multi_prob_gtor_sml_share(const subs_ml_model& smlm)
    : base_t(smlm), _pgoss(smlm) {}

  const prob_gtor_obs_sml_share& get_pgoss() const { return _pgoss; }

private:
  const prob_gtor_obs_sml_share _pgoss;
};



struct root_gtor_sml_share : public multi_prob_gtor_sml_share {
  typedef multi_prob_gtor_sml_share base_t;

  root_gtor_sml_share(const subs_ml_model& smlm) : base_t(smlm) {}

  const prob_t* obs_state_distro_cat(const unsigned cat_no) const;
};



struct bg_gtor_sml_share : public multi_prob_gtor_sml_share {
  typedef multi_prob_gtor_sml_share base_t;

  const notifier& get_obs_info_notifier() const;

  bg_gtor_sml_share(const subs_ml_model& smlm) : base_t(smlm) {}
};


#endif

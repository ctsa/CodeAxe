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

// $Id: lhood_model_prep.h 1182 2008-03-27 01:55:02Z ctsa $

/// \file

#ifndef __LHOOD_MODEL_PREP_H
#define __LHOOD_MODEL_PREP_H

#include "simple_util.h"
#include "subs_ml_types.h"
#include "util/general/uncopyable.h"

#include <memory>
#include <vector>


struct condition_func;
struct group_data;
struct site_data_fastlup;
struct submodel_manager;
struct subs_ml_model;
struct tree_probs;
template <typename T> struct workspace;



/// subset of lhood_model_prep data used in the inner site lhood
/// calculation loop
///
struct lhood_model_prep_cat_site_prob : private uncopyable {

  lhood_model_prep_cat_site_prob(const subs_ml_model& mdl,
                                 const site_data_fastlup& sdf);

  ~lhood_model_prep_cat_site_prob();

  condition_func& cf() { return *_cf; }
  submodel_manager& sm() { return *_sm; }
  tree_probs& tprob() { return *_tprob; }

private:
  std::auto_ptr<condition_func> _cf;
  std::auto_ptr<submodel_manager> _sm;
  std::auto_ptr<tree_probs> _tprob;

public:
  prob_t* site_prob;
  simple_array<workspace<char> > tprob_ws;
  simple_init_matrix<bool> is_adset_using_cat;
  std::vector<std::vector<bool> > cat_site_mask;
};



/// all data that can be initialized and/or malloced independent of
/// the model's parameter state, such that these initializations can be
/// called once before many lhood calls where only the model
/// parameters are allowed to change, (as in during minimization)
///
struct lhood_model_prep : private uncopyable {

  lhood_model_prep(const subs_ml_model& mdl,
                   const site_data_fastlup& sdf);

  ~lhood_model_prep();

  group_data& gd() { return *_gd; }

  const prob_t * const * prob_adset_on_group() const {
    return _prob_adset_on_group.ptr();
  }

private:
  std::auto_ptr<group_data> _gd;

  simple_init_matrix<prob_t> _prob_adset_on_group;

public:
  simple_init_matrix<prob_t> site_cat_mix_prob;
  simple_init_matrix<prob_t> adset_cat_pdistro;
  simple_init_matrix<prob_t> adset_group_cat_pdistro;

  std::vector<smlfloat> ads_group_cat_prior_norm;
  simple_init_matrix<prob_t> group_id_group_cat_pdistro;

  std::vector<prob_t> cat_pdistro;
  std::vector<prob_t> group_cat_pdistro;

  lhood_model_prep_cat_site_prob csp_prep;
};


#endif

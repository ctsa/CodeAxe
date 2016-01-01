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

// $Id: lhood_model_prep_submodel.h 1175 2008-03-26 00:04:19Z ctsa $

/// \file

#ifndef __LHOOD_MODEL_PREP_SUBMODEL_H
#define __LHOOD_MODEL_PREP_SUBMODEL_H

#include "simple_util.h"
#include "site_data_fastlup_core.h"
#include "util/general/uncopyable.h"

#include <memory>
#include <vector>


struct bi_tree;
struct condition_func;
struct site_data_fastlup;
struct subs_ml_model;
struct tree_probs;


struct submodel_data : private uncopyable {

  submodel_data() : is_init(false), n_states(0) {}

  void
  init(const subs_ml_model& mdl,
       const int init_n_states,
       const site_data_fastlup& full_sdf,
       const unsigned* index_translation);

  bool is_init;
  int n_states;
  site_data_fastlup_core sdf;
  std::vector<std::vector<bool> > cat_site_mask;
  simple_init_array<unsigned> full_sdf_index_map;
};




struct submodel_manager : private uncopyable {

  submodel_manager(const subs_ml_model& mdl,
                   const site_data_fastlup& full_sdf);

  condition_func& cf(unsigned i) { return *_cf[i]; }

  submodel_data& submodel(unsigned i) { return _submodel[i]; }
  tree_probs& tprob(unsigned i) { return _tprob[i]; }

  int n_submodels;
  simple_array<simple_init_array<smlfloat> > site_prob;
private:
  simple_array<submodel_data> _submodel;
  simple_array<tree_probs> _tprob;
  simple_array<std::auto_ptr<condition_func> > _cf;
};



#endif

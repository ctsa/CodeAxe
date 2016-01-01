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

// $Id: lhood_model.h 1051 2007-12-05 06:51:26Z ctsa $

/// \file

#ifndef __LHOOD_MODEL_H
#define __LHOOD_MODEL_H

#include "cat_post_prob_mode.h"
#include "simple_util.h"
#include "subs_ml_types.h"


struct lhood_model_prep;
struct site_data_fastlup;
struct site_prob_maker;
struct subs_ml_model;


smlfloat
get_lnp_norm(const subs_ml_model& mdl,
             const site_data_fastlup& sdf);

/// get prob of data given model
///
/// requires a handle to parameter independent data (nlp)
///
void
get_lnprob_from_param_prep(const subs_ml_model& mdl,
                           const site_data_fastlup& full_sdf,
                           lhood_model_prep& nlp,
                           const site_prob_maker& spm,
                           smlfloat& lnp,
                           smlfloat& lnp_norm,
                           const CPPM::index_t cat_post_prob_mode = CPPM::NONE,
                           simple_init_matrix<prob_t>* ppip = 0);

/// get log prob of data given model
///
/// simplest call: this version creates the parameter independent data
/// for you for each call, convenient but less efficient.
///
void
get_lnprob_from_param(const subs_ml_model& mdl,
                      const site_data_fastlup& sdf,
                      smlfloat& lnp,
                      smlfloat& lnp_norm,
                      const CPPM::index_t cat_post_prob_mode = CPPM::NONE,
                      simple_init_matrix<prob_t>* ppip = 0);

#endif

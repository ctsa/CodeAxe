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

// $Id: lhood_model_em.h 1176 2008-03-26 00:17:56Z ctsa $

/// \file

#ifndef __LHOOD_MODEL_EM_H
#define __LHOOD_MODEL_EM_H

#if 0
#include "subs_ml_types.h"

struct bi_tree;
struct estep_dat_t;
struct site_data_fastlup;
struct subs_ml_model;
struct tree_probs;

smlfloat
get_lnp_norm(const subs_ml_model& mdl,
             const estep_dat_t& edat);

// em-version
//
void
get_lnprob_from_param(const subs_ml_model& mdl,
                      const estep_dat_t& edat,
                      smlfloat& lnp,
                      smlfloat& lnp_norm);


/// estep update -- closely related to regular lhood calculation
///
/// !!ALL TPROBS TRANSPOSED!!
void
update_etrans(const bi_tree& tree,
              const tree_probs& tprob,
              const site_data_fastlup& sdf,
              const unsigned n_states,
              estep_dat_t& edat);
#endif
#endif

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

// $Id: lhood_model_em_baz.cc 916 2007-10-12 21:13:20Z ctsa $

/// \file

#include "bi_tree.h"
#include "lhood_model_em_baz.h"
#include "lhood_site_cached_estep_baz.h"
#include "rate_gtor_nscodon_base.h"
#include "site_data_fastlup.h"
#include "subs_ml_model.h"
#include "tree_probs.h"



/// estep update
/// \brief hacky look at the bazykin-kondrashov double nonsynon mut effect
///
void
update_etrans_baz(const subs_ml_model& mdl,
                  const bi_tree& tree,
                  const tree_probs& tprob,
                  const site_data_fastlup& sdf,
                  const unsigned n_states){

#ifdef DEBUG
  // check the distribution and orientation of input probs
  tprob.self_check();
#endif

  const tree_lhood_info ti(tree,tprob,n_states);

  lhood_site_cached_estep_baz(static_cast<const rate_gtor_nscodon_base&>(mdl.get_rate_gtor()),ti,sdf);
}


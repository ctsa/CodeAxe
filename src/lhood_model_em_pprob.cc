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

// $Id: lhood_model_em_pprob.cc 916 2007-10-12 21:13:20Z ctsa $

/// \file

#include "bi_tree.h"
#include "lhood_model_em_pprob.h"
#include "lhood_site_cached_estep_pprob.h"
#include "site_data_fastlup.h"
#include "subs_ml_types.h"
#include "tree_probs.h"



///
void
update_pprob(const bi_tree& tree,
             const tree_probs& tprob,
             const site_data_fastlup& sdf,
             const unsigned n_states,
             pprob_cache& ppc){

#ifdef DEBUG
  // check the distribution and orientation of input probs
  tprob.self_check();
#endif

  const tree_lhood_info ti(tree,tprob,n_states);

  lhood_site_cached_estep_pprob(ti,sdf,ppc);
}

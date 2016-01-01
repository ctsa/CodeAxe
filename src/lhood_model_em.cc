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

// $Id: lhood_model_em.cc 1176 2008-03-26 00:17:56Z ctsa $

/// \file

#if 0
#include "bi_tree.h"
#include "estep_dat_t.h"
#include "lhood_model_em.h"
#include "lhood_site_cached_estep.h"
#include "rate_gtor.h"
#include "site_data_fastlup.h"
#include "subs_ml_model.h"
#ifdef DEBUG
#include "subs_ml_ptol.h"
#endif
#include "subs_ml_types.h"
#include "tree_probs.h"
#include "util/general/log.h"
#include "util/math/array_util.h"
#include "util/math/matrix_util.h"
#include "util/math/test_float.h"

#include <cstring>

#include <ostream>


// tprob values ln'd!!
static
smlfloat
get_mstep_lhood(const tree_probs& tprob,
                const estep_dat_t& edat,
                const unsigned n_states){

  smlfloat mstep_lnp(array_dot(edat.root,tprob.root_prob(),n_states));

  // sum in each branch:
  const unsigned n_branches(edat.n_branches);
  const unsigned n2(n_states*n_states);
  for(unsigned b(0);b<n_branches;++b){
    mstep_lnp += array_dot(edat.trans[b],tprob.branch_prob(b),n2);
  }

  return mstep_lnp;
}




smlfloat
get_lnp_norm(const subs_ml_model& mdl,
             const estep_dat_t& edat){

  smlfloat datasize(0);
  for(unsigned i(0);i<mdl.state_size();++i) datasize += edat.root[i];
  return static_cast<smlfloat>(datasize)*
    std::log(static_cast<smlfloat>(mdl.get_rate_gtor().state_size_conditioned()));
}




// em-version
//
void
get_lnprob_from_param(const subs_ml_model& mdl,
                      const estep_dat_t& edat,
                      smlfloat& lnp,
                      smlfloat& lnp_norm){

  const unsigned n_states(mdl.state_size());
  const unsigned n_branch(mdl.tree().branch_size());

  tree_probs tprob;

  tprob.init_probs(mdl);

  smlfloat adjust_lnp(0.);
#if 0
  // condition probabilities, as specified by the rate_matrix:
  mdl.rate_gtor().cf().mstep_adjust_lnp(mdl.root_gtor().distro(),tprob,adjust_lnp);
#endif

  // min p is setup to handle underflowed transition probabilities at short time
  // intervals
  static const prob_t min_p(1e-100);
  //  static const prob_t log_min_p = std::log(min_p);

  // transpose back to standard form and logify everything in tree_probs
  for(unsigned b(0);b<n_branch;++b){
    matrix_transpose_inplace(tprob.branch_prob(b),n_states);
    array_set_min(tprob.branch_prob(b),min_p,n_states*n_states);
    array_log(tprob.branch_prob(b),n_states*n_states);
  }
  array_set_min(tprob.root_prob(),min_p,n_states);
  array_log(tprob.root_prob(),n_states);

  lnp = get_mstep_lhood(tprob,edat,n_states);

  // condition probabilities, as specified by the rate_matrix:
  lnp += adjust_lnp;

  lnp_norm = lnp / get_lnp_norm(mdl,edat);

  if( lnp_norm <= 2*lnp_norm || is_float_invalid(lnp_norm) ){
    log_os << "FATAL:: Invalid mstep lnp: " << lnp_norm << "\n";
    log_os << "Dumping model:\n";
    mdl.store_state(log_os);
    log_os << "FATAL:: Invalid mstep lnp: " << lnp_norm << "\n";
    abort();
  }
}




/// estep update -- closely related to regular lhood calculation
///
// all tprobs transposed!
void
update_etrans(const bi_tree& tree,
              const tree_probs& tprob,
              const site_data_fastlup& sdf,
              const unsigned n_states,
              estep_dat_t& edat){

  const unsigned n_branches(tree.branch_size());

#ifdef DEBUG
  // check the distribution and orientation of input probs
  tprob.self_check();
#endif

  //zero out eroot/etrans
  memset(edat.root,0,n_states*sizeof(smlfloat));
  for(unsigned i(0);i<n_branches;++i){
    memset(edat.trans[i],0,n_states*n_states*sizeof(smlfloat));
  }

  const tree_lhood_info ti(tree,tprob,n_states);

  lhood_site_cached_estep(ti,sdf,edat);
}
#endif

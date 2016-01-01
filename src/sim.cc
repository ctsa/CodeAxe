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

// $Id: sim.cc 1133 2008-01-29 22:42:35Z ctsa $

/// \file

#include "cat_manager.h"
#include "rate_gtor_nuc.h"
#include "sim.h"
#include "sim_discrete.h"
#include "sim_context.h"
#include "sim_iss.h"

#include "util/general/die.h"



static
bool
test_indy_nuc_model(const subs_ml_model& mdl){

  const rate_gtor_nuc_base* rgn(dynamic_cast<const rate_gtor_nuc_base*>(&mdl.get_rate_gtor()));
  if(rgn==0) return false;

  const cat_manager& cm(mdl.get_cat_manager());
  const unsigned nsmc(cm.typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL));
  for(unsigned i(0);i<nsmc;++i){
    if(rgn->context_model_nuc(i) != CONTEXT_MODEL_NUC::INDY) return false;
  }
  return true;
}



void
simulate_data(const sim_options& sim_opt,
              const subs_ml_model& mdl,
              nuc_seq_data& nsd){

  const RATE_GTOR_MODEL::index_t rgm(mdl.rate_gtor_model());
  const bool is_conditioned_site_model(RATE_GTOR_MODEL::base_overlap_size(rgm) != 0);

  const bool is_indy_nuc_model(test_indy_nuc_model(mdl));

  const bool is_indy_site_model((!is_conditioned_site_model) && is_indy_nuc_model);

  if       ( sim_opt.method == SIM_MODEL::CONTINUOUS){
    // no simpler indy mutation version exists yet:
    simulate_data_context(sim_opt,mdl,nsd);

  } else if( sim_opt.method == SIM_MODEL::DISCRETE ){
    // context simulator is required for conditioned sites to simulate
    // the root distribution, and to handle partial codon selection
    //
    if(is_indy_site_model) simulate_data_discrete(sim_opt,mdl,nsd);
    else                   simulate_data_context(sim_opt,mdl,nsd);

  } else if ( sim_opt.method == SIM_MODEL::ISS ){
    if(! is_indy_site_model)
      pass_away("iss simulation cannot be used with context-dept mutation or conditioned sites");
    simulate_data_indy_site_sampler(sim_opt,mdl,nsd);
  } else {
    die("unknown simulation method");
  }
}

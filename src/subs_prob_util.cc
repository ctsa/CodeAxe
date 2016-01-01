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

// $Id: subs_prob_util.cc 1133 2008-01-29 22:42:35Z ctsa $

/// \file

#include "bi_tree.h"
#include "cat_manager.h"
#include "subs_prob_util.h"
#include "subs_prob_util_core.h"
#include "time_gtor.h"



void
get_all_branch_subs_prob(prob_t** branch_prob,
                         const subs_ml_model& mdl,
                         const rates_func_options_base& bopt,
                         workspace<char> ws[3]) {

  const time_gtor& tg(mdl.get_time_gtor());
  const unsigned n_branches(mdl.tree().branch_size());

  const cat_manager& cm(mdl.get_cat_manager());
  const unsigned n_branch_cat_sets(cm.branch_cat_set_size(bopt.cat));

  simple_array<bool> is_mask(n_branches);

  for(unsigned bcs(0);bcs<n_branch_cat_sets;++bcs){
    const rates_func_options opt(bopt,bcs);
    for(unsigned b(0);b<n_branches;++b){
      is_mask[b]=cm.get_branch_cat_set(bopt.cat,b)!=bcs;
    }
    get_subs_prob(branch_prob,mdl,opt,n_branches,tg.branch_times(opt.cat),ws,is_mask.ptr());
  }
}



void
get_single_branch_subs_prob(prob_t* branch_prob,
                            const smlfloat time,
                            const subs_ml_model& mdl,
                            const rates_func_options_base& bopt,
                            const unsigned branch_id,
                            workspace<char> ws[3]) {

  const unsigned bcs(mdl.get_cat_manager().get_branch_cat_set(bopt.cat,branch_id));
  const rates_func_options opt(bopt,bcs);
  get_subs_prob(&branch_prob,mdl,opt,1,&time,ws);
}

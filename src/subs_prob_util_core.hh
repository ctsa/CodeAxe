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

// $Id: subs_prob_util_core.hh 1061 2007-12-13 22:11:42Z ctsa $

/// \file

#include "branch_time_util.h"
#include "rate_gtor.h"
#include "subs_ml_model_min_options.h"
#include "subs_ml_ptol.h"
#include "subs_prob_util.h"
#include "util/general/die.h"
#include "util/general/log.h"
#include "util/math/matrix_exp_diag.h"
#include "util/math/matrix_exp_series.h"
#include "util/math/prob_util_io.h"


template <typename RandomAccessIterator>
void
get_subs_prob(prob_t** subs_prob,
              const subs_ml_model& mdl,
              const rates_func_options& opt,
              const unsigned n_subs_prob,
              RandomAccessIterator subs_prob_time,
              workspace<char> ws[3],
              const bool* is_mask){

  const unsigned n_states(opt.is_use_submodel ? mdl.submodel_state_size(opt.submodel_no):mdl.state_size());
  const unsigned rate_size(n_states*n_states*sizeof(smlfloat));

  workspace<char>& rate_ws(ws[0]);

  rate_ws.resize(rate_size);

  smlfloat* rates(reinterpret_cast<smlfloat* const>(rate_ws.ptr()));

  mdl.get_rate_gtor().rates(rates,opt);


  if(mdl.opt().is_prefer_diag_expm){
#ifndef USE_LAPACK
    pass_away("Diag expm option cannot be used in non-LAPACK build.");
#else
    // matrix exponential by eigenvector decomposition
    matrix_exp_diag_prepdata pd(n_states);

    matrix_exp_diag_prep(pd,rates);

    for(unsigned b(0);b<n_subs_prob;++b){
      if( is_mask==0 || !is_mask[b]){
        if(is_zero_branch_time(subs_prob_time[b])){
          matrix_identity(subs_prob[b],n_states);
        } else {
          matrix_exp_diag_scale(subs_prob[b],subs_prob_time[b],pd);
        }
      }
    }
#endif
  } else {
    workspace<char>& prep_ws(ws[1]);
    workspace<char>& prep_ws2(ws[2]);
    matrix_exp_series_prepdata<prob_t> mepd(n_states,prep_ws,prep_ws2);

    matrix_exp_series_prep<smlfloat>(mepd,rates);

    for(unsigned b(0);b<n_subs_prob;++b){
      if( is_mask==0 || !is_mask[b]){
        if(is_zero_branch_time(subs_prob_time[b])){
          matrix_identity(subs_prob[b],n_states);
        } else {
          workspace<char>& scale_ws(ws[0]);
          matrix_exp_series_scale<smlfloat>(subs_prob[b],mepd,subs_prob_time[b],scale_ws);
        }
      }
    }
  }

#ifndef NDEBUG
  // check resulting matrices
  for(unsigned b(0);b<n_subs_prob;++b){
    if( is_mask==0 || !is_mask[b]){
      for(unsigned i(0);i<n_states;++i){
        pdistro_check(subs_prob[b]+i*n_states,n_states,SUBS_ML_PTOL);
      }
    }
  }
#endif
}

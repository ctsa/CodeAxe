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

// $Id: subs_ml_fisher.cc 1111 2008-01-26 01:11:52Z ctsa $

/// \file

#include "lhood_model.h"
#include "subs_ml_fisher.h"
#include "util/general/log.h"
#include "util/math/deriv_util.h"

#include <ostream>


struct subs_ml_param_func {

  subs_ml_param_func(const site_data_fastlup& init_sdf,
                     const subs_ml_model& init_mdl,
                     const unsigned init_param_no)
    : sdf(init_sdf), mdl(init_mdl), param_no(init_param_no) {}

  smlfloat operator()(const smlfloat val) const {

    const unsigned n(mdl.param_size());

    simple_array<bool> is_train(n);
    mdl.is_train_param_state(is_train.begin());

    simple_array<smlfloat> pcopy(n);
    mdl.param_state(pcopy.begin());

    subs_ml_model mdl_copy(mdl);

    if(is_train[param_no]) {
      is_train[param_no] = false;
      mdl_copy.set_is_train_param_state(is_train.begin());
    }

    pcopy[param_no] = val;
    mdl_copy.set_param_state(pcopy.begin());

    smlfloat lnp,lnp_norm;
    get_lnprob_from_param(mdl_copy,sdf,lnp,lnp_norm);

    return lnp;
  }

  const site_data_fastlup& sdf;
  const subs_ml_model& mdl;
  const unsigned param_no;
};




struct subs_ml_2param_func {

  subs_ml_2param_func(const site_data_fastlup& init_sdf,
                      const subs_ml_model& init_mdl,
                      const unsigned init_param_no1,
                      const unsigned init_param_no2)
    : sdf(init_sdf), mdl(init_mdl),
      param_no1(init_param_no1), param_no2(init_param_no2) {}

  smlfloat operator()(const smlfloat val1,
                      const smlfloat val2) const {

    const unsigned n(mdl.param_size());

    simple_array<bool> is_train(n);
    mdl.is_train_param_state(is_train.begin());

    simple_array<smlfloat> pcopy(n);
    mdl.param_state(pcopy.begin());

    subs_ml_model mdl_copy(mdl);

    if(is_train[param_no1] || is_train[param_no2]){
      is_train[param_no1] = false;
      is_train[param_no2] = false;
      mdl_copy.set_is_train_param_state(is_train.begin());
    }

    pcopy[param_no1] = val1;
    pcopy[param_no2] = val2;
    mdl_copy.set_param_state(pcopy.begin());

    smlfloat lnp,lnp_norm;
    get_lnprob_from_param(mdl_copy,sdf,lnp,lnp_norm);

    return lnp;
  }

  const site_data_fastlup& sdf;
  const subs_ml_model& mdl;
  const unsigned param_no1;
  const unsigned param_no2;
};




void
get_fisher_matrix(const site_data_fastlup& sdf,
                  const subs_ml_model& mdl,
                  smlfloat* fmat) {

  static const smlfloat delta(1e-4);

  const unsigned n(mdl.param_size());
  const unsigned tn(mdl.param_size(PARAM_VIEW::TRAINABLE));

  simple_array<smlfloat> pcopy(n);
  mdl.param_state(pcopy.begin());

  simple_array<bool> is_train(n);
  mdl.is_train_param_state(is_train.begin());

  smlfloat lnp,lnp_norm;
  get_lnprob_from_param(mdl,sdf,lnp,lnp_norm);

  for(unsigned i(0),ti(0),tij(0);i<n;++i){
    if(! is_train[i]) continue;
    for(unsigned j(i),tj(ti);j<n;++j){
      if(! is_train[j]) continue;

      if((tij+1)%50==0) log_os << ".";
      if((tij+1)%3000==0) log_os << "\n";
      if(i==j){
        subs_ml_param_func f(sdf,mdl,j);
        fmat[tj+tn*ti] = -centered_2nd_deriv_estimate(pcopy[j],lnp,f,delta);
      } else {
        subs_ml_2param_func f_j(sdf,mdl,i,j);
        fmat[tj+tn*ti] = -centered_2partial_deriv_estimate(pcopy[i],pcopy[j],f_j,delta);
        fmat[ti+tn*tj] =  fmat[tj+tn*ti];
      }
      tj++; tij++;
    }
    ti++;
  }
  log_os << "\n";
}

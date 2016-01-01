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

// $Id: bg_gtor.cc 1133 2008-01-29 22:42:35Z ctsa $

/// \file

#include "bg_gtor.h"
#include "cat_manager.h"
#include "subs_ml_print_util.h"
#include "util/general/die.h"
#include "util/general/io_util.h"
#include "util/general/log.h"
#include "util/math/random_util.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>



bg_gtor::
bg_gtor(const self_t& s,
        const bg_gtor_sml_share& smls)
  : base_t(s), _smls(smls), _ptor(s._ptor.size()), _data(s._data) {
  const unsigned ds(distro_cat_size());
  if       (ptype()==PT_PARAM){
    for(unsigned i(0);i<ds;++i) _ptor[i]=s.distro_cat_ptorp(i).clone();
  } else if(ptype()==PT_OBS){
    for(unsigned i(0);i<ds;++i) _ptor[i]=s.distro_cat_ptoro(i).clone(_smls.get_pgoss());
  } else {
    die("unknown ptype");
  }
  register_param();
  observe_notifier(_smls.get_obs_info_notifier());
}



void
bg_gtor::
reset_param(const PARAM_INIT_TYPE::index_t pinit){
  if     (ptype()==PT_OBS) return;
  else if(ptype()!=PT_PARAM) {
    die("unknown ptype");
  }

  const unsigned ds(distro_cat_size());
  for(unsigned i(0);i<ds;++i){
    ptorp_t& p(distro_cat_ptorp(i));

    if       ( pinit == PARAM_INIT_TYPE::RANDOM ){
      static const smlfloat unorm_min(1.e-1);
      const unsigned ps(p.param_size());
      simple_array<smlfloat> param(ps);
      for(unsigned j(0);j<ps;++j){ param[j] = random_uniform()+unorm_min; }
      p.set_param_state(param.begin());
    } else if( pinit == PARAM_INIT_TYPE::START ){
      const unsigned ss(p.state_size());
      simple_array<smlfloat> distro(ss);
      pdistro_unif(distro.begin(),distro.end());
      p.set_pdistro(distro.begin());
    } else { die("unknown param_init_type"); }
  }
  notify_observers(EVENT_TYPE::PDISTRO_CHANGE);
}



const char* const section_id = "bg_gtor";
const char* const end_label = "END";

const iliner il(section_id);



void
bg_gtor::
store_state(std::ostream& os) const {

  os << section_id << " subs_rate_bg_model " << static_cast<int>(_data.bgo.bg) << "\n";

  if(ptype()!=PT_OBS){
    const unsigned ds(distro_cat_size());
    for(unsigned i(0);i<ds;++i){
      std::ostringstream oss;
      oss << section_id << " bg_cat_" << i;
      distro_cat_ptorp(i).store_state(oss.str().c_str(),os);
    }
  }
  os << section_id << " " << end_label << "\n";
}



void
bg_gtor::
load_state(std::istream& is) {

  il.advance(is,"subs_rate_bg_model");
  _data.bgo.bg=static_cast<SUBS_RATE_BG_MODEL::index_t>(read_int(is));

  init();

  if(ptype()!=PT_OBS){
    const unsigned ds(distro_cat_size());
    for(unsigned i(0);i<ds;++i){
      ptorp_t& p(distro_cat_ptorp(i));
      const unsigned ps(p.param_size());
      simple_array<smlfloat> param(ps);
      simple_array<bool> itp(ps);
      for(unsigned j(0);j<ps;++j){
        il.advance(is);
        bool btmp;
        is >> btmp >> param[j];
        itp[j] = btmp;
      }
      p.set_param_state(param.begin());
      p.set_is_train_param_state(itp.begin());
    }
    notify_observers(EVENT_TYPE::PDISTRO_CHANGE);
  }
  il.advance(is,end_label);
}




// determine the minimum number of distro cats required to represent
// the model category structure if the bg is dependent on observed
// state distributions. For parameterized bg models it's much
// simpler: n_distro_cats = n_bg_cats.
//
void
bg_gtor::
setup_distro_cats(){

  const unsigned n_cats(cm().cat_size());
  _data.distro_cat_map.resize(n_cats);

  if(ptype()==PT_PARAM){
    for(unsigned i(0);i<n_cats;++i){
      _data.distro_cat_map[i]=cm().typed_cat_no_from_cat_no(i,CAT_PARAM_TYPE::SUBS_RATE_BG);
    }
    _data.distro_cat_size=cm().typed_cat_size(CAT_PARAM_TYPE::SUBS_RATE_BG);

  } else if(ptype()==PT_OBS){
    using namespace SUBS_RATE_BG_MODEL;

    for(unsigned i(0);i<n_cats;++i){
      _data.distro_cat_map[i]=cm().typed_cat_no_from_cat_no(i,CAT_PARAM_TYPE::OBS);
    }
    _data.distro_cat_size=cm().typed_cat_size(CAT_PARAM_TYPE::OBS);

  } else {
    die("Unknown ptype");
  }
}



#include "prob_gtor_param_state_conversion.h"
#include "prob_gtor_param_simple.h"
#include "prob_gtor_obs_derived.h"



bg_gtor::ptorp_t*
bg_gtor::
prob_gtor_param_factory() const {

  const unsigned n_states(state_size());
  const RATE_GTOR_MODEL::index_t ragm(_smls.rate_gtor_model());

  const bool is_c4_site_model(ragm == RATE_GTOR_MODEL::C4PRE ||
                              ragm == RATE_GTOR_MODEL::C4POST);

  const bool is_c4_pre(ragm==RATE_GTOR_MODEL::C4PRE);


  if       (_data.bgo.bg == SUBS_RATE_BG_MODEL::CODON_DINUC){
    if       (is_c4_site_model){
      return new prob_gtor_param_nsc4_seq_stat_nx1n(is_c4_pre);
    } else {
      pass_away("bg codon-dinuc model only valid with c4-pre/c4-post site models");
    }
  } else if(_data.bgo.bg == SUBS_RATE_BG_MODEL::FULL){
    if(is_c4_site_model){
      return new prob_gtor_param_nsc4_seq_stat(is_c4_pre);
    } else {
      return new prob_gtor_param_simple(n_states);
    }
  } else {
    pass_away("unrecognized bg-model type");
  }
}



bg_gtor::ptoro_t*
bg_gtor::
prob_gtor_obs_factory(const unsigned distro_cat_no) const {

  using namespace SUBS_RATE_BG_MODEL;

  unsigned cat_no(0);
  const unsigned n_cats(cm().cat_size());
  for(unsigned c(0);c<n_cats;++c){
    if(_data.distro_cat_map[c] == distro_cat_no){
      cat_no=c;
      break;
    }
    if((c+1)==n_cats){
      die("unmapped distro_cat_no");
    }
  }

  if       (_data.bgo.bg == OBS_AVG){
    return new prob_gtor_obs_state_avg(_smls.get_pgoss(),cat_no);
  } else {
    die("unkown bg_model in prob_gtor_obs_factory");
  }
}

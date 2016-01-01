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

// $Id: obs_info.cc 1121 2008-01-28 22:32:57Z ctsa $

/// \file

#include "obs_info.h"
#include "subs_ml_print_util.h"
#include "util/general/die.h"
#include "util/general/io_util.h"
#include "util/general/log.h"
#include "util/math/prob_util.h"
#include "util/math/random_util.h"

#include <algorithm>
#include <iomanip>
#include <iostream>



void
obs_info::
update_distros(){

  if(_source==NONE) die("Invalid obs_info state.");
  const unsigned n_obs_cats(_data.cat_seq_state_source.dim1());
  const unsigned n_seqs(_data.cat_seq_state_source.dim2());
  const unsigned n_states(_data.cat_seq_state_source.dim3());

  {
    _data.state_distro.set_val(0.);
    prob_t* distro(_data.state_distro.ptr());

    for(unsigned j(0);j<n_seqs;++j){
      prob_t* sdistro(_data.seq_state_distro[j]);
      for(unsigned k(0);k<n_states;++k){
        sdistro[k]=0.;
        for(unsigned i(0);i<n_obs_cats;++i){
          const prob_t val(_data.cat_seq_state_source[i][j][k]);
          sdistro[k] += val;
          distro[k] += val;
        }
      }
      pdistro_norm_safe(sdistro,sdistro+n_states);
    }
    pdistro_norm_safe(distro,distro+n_states);
  }

  for(unsigned i(0);i<n_obs_cats;++i){
    prob_t* cdistro(_data.state_distro_obs_cat[i]);
    std::fill(cdistro,cdistro+n_states,0.);
    for(unsigned j(0);j<n_seqs;++j){
      prob_t* scdistro(_data.seq_state_distro_obs_cat[i][j]);
      for(unsigned k(0);k<n_states;++k){
        const prob_t val(_data.cat_seq_state_source[i][j][k]);
        scdistro[k] = val;
        cdistro[k] += val;
      }
      pdistro_norm_safe(scdistro,scdistro+n_states);
    }
    pdistro_norm_safe(cdistro,cdistro+n_states);
  }

  notify_observers(EVENT_TYPE::PDISTRO_CHANGE);
}



template <typename T>
void
obs_info::
set_obs_cat_seq_state_int(T const * const * const * c){
  const unsigned n_obs_cats(_data.cat_seq_state_source.dim1());
  const unsigned n_seqs(_data.cat_seq_state_source.dim2());
  const unsigned n_states(_data.cat_seq_state_source.dim3());

  for(unsigned i(0);i<n_obs_cats;++i){
    for(unsigned j(0);j<n_seqs;++j){
      for(unsigned k(0);k<n_states;++k){
        _data.cat_seq_state_source[i][j][k] = static_cast<smlfloat>(c[i][j][k]);
      }
    }
  }

  update_distros();
}



void
obs_info::
set_obs_cat_seq_state_distro(prob_t const * const * const * c){
  _source=DISTRO;
#ifdef DEBUG
  /// \todo add distro check for this case:
  //  pdistro_check
#endif
  set_obs_cat_seq_state_int(c);
}

void
obs_info::
set_obs_cat_seq_state_counts(smlfloat const * const * const * c){
  _source=COUNTS;
  set_obs_cat_seq_state_int(c);
}



const char* const section_id = "obs_info";
const char end_label[] = "END";
const char* const oss_tag = "obs_state_source";

const iliner il(section_id);



void
obs_info::
store_state(std::ostream& os) const {

  const unsigned n_obs_cats(_data.cat_seq_state_source.dim1());
  const unsigned n_seqs(_data.cat_seq_state_source.dim2());
  const unsigned n_states(_data.cat_seq_state_source.dim3());

  os << section_id << " " << oss_tag << " " << _source << "\n";
  os << std::setprecision(SUBS_ML_PRINT_MODEL_PRECISION);
  for(unsigned k(0);k<n_obs_cats;++k){
    for(unsigned i(0);i<n_states;++i){
      os << section_id << " obs_" << k << "_seq_state_" << i;
      for(unsigned j(0);j<n_seqs;++j) os << " " << _data.cat_seq_state_source[k][j][i];
      os << "\n";
    }
  }
  os << section_id << " " << end_label << "\n";
}



void
obs_info::
reset_param(const PARAM_INIT_TYPE::index_t pinit){
  _source=DISTRO;
  if(pinit == PARAM_INIT_TYPE::RANDOM){
    static const smlfloat unorm_min(1.e-1);
    const unsigned d1(_data.cat_seq_state_source.dim1());
    const unsigned d2(_data.cat_seq_state_source.dim2());
    const unsigned d3(_data.cat_seq_state_source.dim3());
    for(unsigned i(0);i<d1;++i){
      for(unsigned j(0);j<d2;++j){
        smlfloat* p(_data.cat_seq_state_source[i][j]);
        for(unsigned k(0);k<d3;++k) p[k] = random_uniform()+unorm_min;
        pdistro_norm(p,p+d3);
      }
    }
  } else if(pinit == PARAM_INIT_TYPE::START){
    _data.cat_seq_state_source.set_val(0);
  } else { die("unknown param_init_type"); }

  update_distros();
}



void
obs_info::
load_state(std::istream& is) {

  const unsigned n_obs_cats(_data.cat_seq_state_source.dim1());
  const unsigned n_seqs(_data.cat_seq_state_source.dim2());
  const unsigned n_states(_data.cat_seq_state_source.dim3());

  il.advance(is,oss_tag);
  _source=static_cast<source_type>(read_int(is));

  for(unsigned k(0);k<n_obs_cats;++k){
    for(unsigned i(0);i<n_states;++i){
      il.advance(is);
      for(unsigned j(0);j<n_seqs;++j) is >> _data.cat_seq_state_source[k][j][i];
    }
  }
  update_distros();
  il.advance(is,end_label);
}

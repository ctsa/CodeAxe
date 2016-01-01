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

// $Id: root_gtor.cc 1118 2008-01-28 20:32:23Z ctsa $

/// \file

#include "cat_manager.h"
#include "prob_gtor_param_gc_shifter.h"
#include "root_gtor.h"
#include "subs_ml_print_util.h"
#include "util/general/die.h"
#include "util/general/io_util.h"
#include "util/general/log.h"
#include "util/math/random_util.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>



root_gtor::
root_gtor(const self_t& s,
          const root_gtor_sml_share& smls)
  : base_t(s), _smls(smls), _ptor(s._ptor.size()), _data(s._data) {
  const unsigned ds(distro_cat_size());
  if       (ptype()==PT_PARAM){
    for(unsigned i(0);i<ds;++i) _ptor[i]=s.distro_cat_ptorp(i).clone();
  } else if(ptype()==PT_PARAM_LINKED){
    for(unsigned i(0);i<ds;++i){
      if(i==0){
        _ptor[i]=s.distro_cat_ptorp(i).clone();
      } else {
        _ptor[i]=new prob_gtor_param_gc_shifter(static_cast<const prob_gtor_param_gc_shifter&>(s.distro_cat_ptorp(i)),
                                                distro_cat_ptorp(0));
      }
    }
  } else {
    for(unsigned i(0);i<ds;++i) _ptor[i]=s.distro_cat_ptoro(i).clone(_smls.get_pgoss());
  }
  register_param();
}



void
root_gtor::
reset_param(const PARAM_INIT_TYPE::index_t pinit){
  if(ptype()==PT_OBS) return;

  const unsigned ds(distro_cat_size());
  for(unsigned i(0);i<ds;++i){
    ptorp_t& p(distro_cat_ptorp(i));

    if(i==0 || ptype()==PT_PARAM){
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
    } else if(ptype()==PT_PARAM_LINKED){
      if       ( pinit == PARAM_INIT_TYPE::RANDOM ){
        const smlfloat ran_range(10.);
        const smlfloat param((random_uniform()-.5)*ran_range);
        p.set_param_state(&param);
      } else if( pinit == PARAM_INIT_TYPE::START ){
        const smlfloat param(0.);
        p.set_param_state(&param);
      } else { die("unknown param_init_type"); }
    } else {
      die("unknown ptype");
    }
  }
}



const char* const section_id = "root_gtor";
const char* const end_label = "END";

const iliner il(section_id);



void
root_gtor::
store_state(std::ostream& os) const {

  os << section_id << " root_gtor_model " << static_cast<int>(_data.rogm) << "\n";

  if(ptype()!=PT_OBS){
    const unsigned ds(distro_cat_size());
    for(unsigned i(0);i<ds;++i){
      std::ostringstream oss;
      oss << section_id << " root_cat_" << i;
      distro_cat_ptorp(i).store_state(oss.str().c_str(),os);
    }
  }
  os << section_id << " " << end_label << "\n";
}



void
root_gtor::
load_state(std::istream& is) {

  il.advance(is,"root_gtor_model");
  _data.rogm=static_cast<ROOT_GTOR_MODEL::index_t>(read_int(is));

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
  }
  il.advance(is,end_label);
}



// determine the minimum number of distro cats required to represent
// the model category structure if the root is dependent on model
// parameters/observed state distributions. For parameterized root
// models it's much simpler: n_distro_cats = n_root_cats.
//
void
root_gtor::
setup_distro_cats(){

  const unsigned n_cats(cm().cat_size());
  const unsigned n_branches(cm().branch_size());
  _data.distro_cat_map.resize(n_cats);

  if(ptype()==PT_PARAM || ptype()==PT_PARAM_LINKED){
    for(unsigned i(0);i<n_cats;++i){
      _data.distro_cat_map[i]=cm().typed_cat_no_from_cat_no(i,CAT_PARAM_TYPE::ROOT);
    }
    _data.distro_cat_size=cm().typed_cat_size(CAT_PARAM_TYPE::ROOT);

  } else if(ptype()==PT_OBS){
    using namespace ROOT_GTOR_MODEL;

    const index_t r(_data.rogm);

    if       (r==OBS_AVG){
      for(unsigned i(0);i<n_cats;++i){
        _data.distro_cat_map[i]=cm().typed_cat_no_from_cat_no(i,CAT_PARAM_TYPE::OBS);
      }
      _data.distro_cat_size=cm().typed_cat_size(CAT_PARAM_TYPE::OBS);

    } else if(r==OBS_TIME_AVG || r==LSPROB){
      std::map<std::vector<unsigned>,unsigned> cat_reduce;
      std::vector<unsigned> xcats;

      _data.distro_cat_size=0;
      for(unsigned c(0);c<n_cats;++c){
        xcats.clear();

        if(r==OBS_TIME_AVG){
          xcats.push_back(cm().typed_cat_no_from_cat_no(c,CAT_PARAM_TYPE::OBS));
          xcats.push_back(cm().typed_cat_no_from_cat_no(c,CAT_PARAM_TYPE::TIME));
        } else {
          for(unsigned pt(0);pt<CAT_PARAM_TYPE::SIZE;++pt){
            if(pt==CAT_PARAM_TYPE::ROOT) continue;
            for(unsigned mt(0);mt<CAT_MIX_TYPE::SIZE;++mt){
              for(unsigned b(0);b<n_branches;++b){
                xcats.push_back(cm().typed_cat_no_from_cat_no_y_branch_id(c,b,
                                                                          static_cast<const CAT_PARAM_TYPE::index_t>(pt),
                                                                          static_cast<const CAT_MIX_TYPE::index_t>(mt)));
              }
            }
          }
        }

        if( cat_reduce.find(xcats) == cat_reduce.end() ){
          cat_reduce[xcats] = _data.distro_cat_size++;
        }
        _data.distro_cat_map[c] = cat_reduce[xcats];
      }
    } else {
      die("invalid root state");
    }
  } else {
    die("Unknown ptype");
  }
}



#include "prob_gtor_param_state_conversion.h"
#include "prob_gtor_param_simple.h"
#include "prob_gtor_obs_derived.h"



root_gtor::ptorp_t*
root_gtor::
prob_gtor_param_factory(const unsigned distro_cat_no) const {

  const unsigned n_states(state_size());
  const RATE_GTOR_MODEL::index_t ragm(_smls.rate_gtor_model());

  const bool is_conditioned_site_model(RATE_GTOR_MODEL::base_overlap_size(ragm) != 0);

  const bool is_c4_site_model(ragm == RATE_GTOR_MODEL::C4PRE ||
                              ragm == RATE_GTOR_MODEL::C4POST);

  const bool is_c4_pre(ragm==RATE_GTOR_MODEL::C4PRE);


  if       (_data.rogm == ROOT_GTOR_MODEL::CODON){
    if       (is_c4_site_model){
      return new prob_gtor_param_nsc4_from_nscodon(is_c4_pre);
    } else if(ragm==RATE_GTOR_MODEL::C5){
      return new prob_gtor_param_nsc5_from_nscodon();
    } else {
      pass_away("root codon model only valid with c4-pre/c4-post/c5 site models");
    }
  } else if(_data.rogm == ROOT_GTOR_MODEL::CODON_DINUC){
    if       (is_c4_site_model){
      return new prob_gtor_param_nsc4_seq_stat_nx1n(is_c4_pre);
    } else {
      pass_away("root codon-dinuc model only valid with c4-pre/c4-post site models");
    }
  } else if(_data.rogm == ROOT_GTOR_MODEL::CODON_TRINUC){
    if       (is_c4_site_model){
      return new prob_gtor_param_nsc4_seq_stat_nx2n(is_c4_pre);
    } else {
      pass_away("root codon-trinuc model only valid with c4-pre/c4-post site models");
    }
  } else if(_data.rogm == ROOT_GTOR_MODEL::NUC){
    if       (ragm==RATE_GTOR_MODEL::DINUC){
      return new prob_gtor_param_dinuc_from_nuc();
    } else if(ragm==RATE_GTOR_MODEL::TRINUC){
      return new prob_gtor_param_trinuc_from_nuc();
    } else if(ragm==RATE_GTOR_MODEL::CODON){
      return new prob_gtor_param_nscodon_from_nuc();
    } else if(is_c4_site_model){
      return new prob_gtor_param_nsc4_from_nuc();
    } else {
      pass_away("invalid site model for root nuc model");
    }
  } else if(_data.rogm == ROOT_GTOR_MODEL::NUC_POS){
    if(ragm==RATE_GTOR_MODEL::CODON){
      return new prob_gtor_param_nscodon_from_nuc_pos();
    } else if(is_c4_site_model){
      return new prob_gtor_param_nsc4_from_nuc_pos(is_c4_pre);
    } else {
      pass_away("invalid site model root nuc-pos model");
    }
  } else if(_data.rogm == ROOT_GTOR_MODEL::FULL_GC_SHIFT){
    if(is_conditioned_site_model) pass_away("gc_shift root requires non-overlapping sites");

    if(distro_cat_no == 0){
      return new prob_gtor_param_simple(n_states);
    } else {
      return new prob_gtor_param_gc_shifter(distro_cat_ptorp(0),
                                            RATE_GTOR_MODEL::convert_to_site_model(ragm));
    }
  } else if(_data.rogm == ROOT_GTOR_MODEL::FULL){
    if(is_c4_site_model){
      return new prob_gtor_param_nsc4_seq_stat(is_c4_pre);
    } else {
      return new prob_gtor_param_simple(n_states);
    }
  } else {
    pass_away("unrecognized root-model type");
  }
}



root_gtor::ptoro_t*
root_gtor::
prob_gtor_obs_factory(const unsigned distro_cat_no) const {

  using namespace ROOT_GTOR_MODEL;

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

  const prob_gtor_obs_sml_share& pgoss(_smls.get_pgoss());

  if       (_data.rogm == OBS_AVG){
    return new prob_gtor_obs_state_avg(pgoss,cat_no);
  } else if(_data.rogm == OBS_TIME_AVG){
    return new prob_gtor_obs_state_time_avg(pgoss,cat_no);
  } else if(_data.rogm == LSPROB){
    return new prob_gtor_obs_root_lsprob(pgoss,cat_no);
  } else {
    die("unkown root_gtor_model in prob_gtor_obs_factory");
  }
}

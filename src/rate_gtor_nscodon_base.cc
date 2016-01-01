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

// $Id: rate_gtor_nscodon_base.cc 1183 2008-03-27 02:12:11Z ctsa $

/// \file

#include "cat_manager.h"
#include "nsaa_maps.h"
#include "rate_gtor_nscodon_base.h"
#include "rate_gtor_nscodon_base_util.h"
#include "rate_gtor_util.h"
#include "substk_exception.h"
#include "util/bio/bioseq_util_pdf_conversions.h"
#include "util/general/io_util.h"
#include "util/general/log.h"
#include "util/math/random_util.h"

#include <cassert>
#include <cstdio> // for snprintf

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>



// instantiate static class data
const rate_gtor_nscodon_base::aa_exchange_param_info rate_gtor_nscodon_base::_aaex;
const rate_gtor_nscodon_base::aa_codon_map_info rate_gtor_nscodon_base::_aacodon;



rate_gtor_nscodon_base::
aa_exchange_param_info::
aa_exchange_param_info(){
  make_reachable_aa_matrix(is_1nuc_aa_exchange);

  // build aa parameter offset maps based on single nuc exchanges:
  //
  symm_size = 0;
  asymm_size = 0;

  {
    static const unsigned AA2(NSAA::SIZE*NSAA::SIZE);
    std::fill(asymm_aatoparam_map,asymm_aatoparam_map+AA2,0);
    std::fill(symm_aatoparam_map,symm_aatoparam_map+AA2,0);
  }
  asymm_paramtoaa_map.clear();
  symm_paramtoaa_map.clear();

  for(unsigned i(0);i<NSAA::SIZE;++i){
    for(unsigned j(0);j<NSAA::SIZE;++j){
      const unsigned aaindex(j+i*NSAA::SIZE);
      if(! is_1nuc_aa_exchange[aaindex] ) continue;

      asymm_aatoparam_map[aaindex] = asymm_size;
      asymm_paramtoaa_map.push_back(aaindex);
      asymm_size += 1;
    }
    for(unsigned j(i);j<NSAA::SIZE;++j){
      const unsigned aaindex(j+i*NSAA::SIZE);
      if(! is_1nuc_aa_exchange[aaindex] ) continue;

      symm_aatoparam_map[aaindex] = symm_size;
      symm_aatoparam_map[i+j*NSAA::SIZE] = symm_size;
      symm_paramtoaa_map.push_back(aaindex);
      symm_size += 1;
    }
  }
}


rate_gtor_nscodon_base::
aa_codon_map_info::
aa_codon_map_info(){
  // find all codons for each amino acid & create a codon order remap
  // which causes the translated amino acids to fall into continuous
  // blocks in the amino acid index order
  //
  std::vector<std::vector<NSCODON::index_t> > aa_codons(NSAA::SIZE);
  for(unsigned i(0);i<NSCODON::SIZE;++i){
    const NSCODON::index_t c(static_cast<NSCODON::index_t>(i));
    const NSAA::index_t aa(codon_trans_known(c));
    aa_codons[aa].push_back(c);
  }

  unsigned ci(0);
  for(unsigned i(0);i<NSAA::SIZE;++i){
    const unsigned aas(aa_codons[i].size());
    aa_codon_size[i]=aas;
    aa_codon_start_offset[i]=ci;
    for(unsigned j(0);j<aas;++j){
      nscodon_aa2normal_order_map[ci] = aa_codons[i][j];
      ci++;
    }
  }

  for(unsigned i(0);i<NSCODON::SIZE;++i){
    nscodon_normal2aa_order_map[nscodon_aa2normal_order_map[i]] = i;
  }
}


static
bool
is_sm_scaled(const SELECT_MODEL::index_t sm){
  using namespace SELECT_MODEL;
  switch(sm){
  case NONE:
  case SINGLE: return false;
  default:     return true;
  }
}


unsigned
rate_gtor_nscodon_base::
paramsize_select_model_unscaled(const SELECT_MODEL::index_t sm) {

  using namespace SELECT_MODEL;

  switch(sm){
  case NONE:         return 0;
  case SINGLE:       return 1;
  case HP:           return 4;
  case FROM:
  case TO:           return NSAA::SIZE;
  case FROM_TO:      return NSAA::SIZE*2;
  case SYMM:         return _aaex.symm_size;
  case FROM_SYMM:
  case SYMM_TO:      return _aaex.symm_size + NSAA::SIZE;
  case FROM_SYMM_TO: return _aaex.symm_size + NSAA::SIZE*2;
  case ASYMM:        return _aaex.asymm_size;
  default:     die("Unknown selection model");
  }
}


unsigned
rate_gtor_nscodon_base::
paramsize_select_model(const SELECT_MODEL::index_t sm) {

  const unsigned scale(is_sm_scaled(sm) ? 1 : 0);
  return paramsize_select_model_unscaled(sm)+scale;
}


void
rate_gtor_nscodon_base::
init_class_storage(){

  _data.copt.set_n_select_matrix_cats(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_MATRIX));

  resize(param_end());

  {
    const unsigned st(base_t::param_end());
    const unsigned sts(param_end());
    for(unsigned i(st);i<sts;++i){
      _is_free_block_start[i] = false;
      _is_free_block_stop[i] = false;
    }
  }

  if(_data.copt.codon_bias_model == CODON_BIAS_MODEL::SYNON_RATIO ||
     _data.copt.codon_bias_model == CODON_BIAS_MODEL::SYNON_NORM_RATIO ||
     _data.copt.codon_bias_model == CODON_BIAS_MODEL::SYNON_NORM_ABS){
    const unsigned pst(paramset_start(PARAMSET_CODON_BIAS));
    for (unsigned i(0);i<NSAA::SIZE;++i){
      const unsigned st(pst+_aacodon.aa_codon_start_offset[i]);
      const unsigned s(_aacodon.aa_codon_size[i]);
      if(s){
        _is_free_block_start[st] = true;
        _is_free_block_stop[st+s-1] = true;
      }
    }
  }

  /// \todo these paramset blocks are probably not accurate in the new cat scheme:
  set_paramset_dependent_block(PARAMSET_SITE_SELECT_CATS);
  set_paramset_dependent_block(PARAMSET_GROUP_SELECT_CATS);
}



void
rate_gtor_nscodon_base::
init_param_locks(const rate_gtor_nscodon_lock_options& lopt){

  // lock parameters which won't be trained
  {
    CAT_PARAM_TYPE::index_t pt(CAT_PARAM_TYPE::SEL_STRENGTH);
    if(lopt.is_lock_site_select_cats){
      set_paramset_train(PARAMSET_SITE_SELECT_CATS,false);
    } else {
      CAT_MIX_TYPE::index_t mt(CAT_MIX_TYPE::SITE);
      lock_1param_norm_sets(PARAMSET_SITE_SELECT_CATS,pt,mt,
                            cm().sel_strength_norm_set_size(mt),
                            cm().sel_strength_norm_set_from_sel_strength_cat(mt));
    }

    if(lopt.is_lock_group_select_cats){
      set_paramset_train(PARAMSET_GROUP_SELECT_CATS,false);
    } else {
      CAT_MIX_TYPE::index_t mt(CAT_MIX_TYPE::GROUP);
      lock_1param_norm_sets(PARAMSET_GROUP_SELECT_CATS,pt,mt,
                            cm().sel_strength_norm_set_size(mt),
                            cm().sel_strength_norm_set_from_sel_strength_cat(mt));
    }
  }

  const unsigned nsmc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_MATRIX));
  for(unsigned i(0);i<nsmc;++i){
    using namespace SELECT_MODEL;
    const index_t s(select_model(i));
    if(! (s==NONE || s==SINGLE)){
      _is_train_param[paramset_aa_select_matrix_cat_start(i)] = false;
    }
  }

  if(_data.copt.codon_bias_model == CODON_BIAS_MODEL::SYNON_RATIO ||
     _data.copt.codon_bias_model == CODON_BIAS_MODEL::SYNON_NORM_RATIO){
    const unsigned st(paramset_start(PARAMSET_CODON_BIAS));
    for(unsigned i(0);i<NSAA::SIZE;++i){
      if(_aacodon.aa_codon_size[i]==1){
        _is_train_param[st+_aacodon.aa_codon_start_offset[i]] = false;
      }
    }
  }
}




void
rate_gtor_nscodon_base::
rates(smlfloat rts[],
      const rates_func_options& opt) const {
  rates_nscodon_context(*this,rts,opt);
}


void
rate_gtor_nscodon_base::
rates_variance(irv_t<smlfloat> rts[],
               const rates_func_options& opt) const {
  rates_nscodon_context(*this,rts,opt);
}


void
rate_gtor_nscodon_base::
state_pdistro_to_nscodon_pdistro(const prob_t* state_pdistro,
                                 prob_t* nscodon_pdistro) const {
  const BIO_PDISTRO::index_t i(BIO_PDISTRO::convert_from_site_model(site_model()));
  bio_pdistro_convert(i,BIO_PDISTRO::NSCODON,state_pdistro,nscodon_pdistro);
}


static
bool
is_cbias_pseudo_param(const CODON_BIAS_MODEL::index_t& cbm){
  using namespace CODON_BIAS_MODEL;
  switch(cbm){
  case SYNON_RATIO:
  case SYNON_NORM_RATIO: return true;
  default:               return false;
  }
}



unsigned
rate_gtor_nscodon_base::
pseudo_param_size() const {
  unsigned n_pps(0);

  // pseudo-parameter to power-scale all codon bias parameters
  if(is_cbias_pseudo_param(_data.copt.codon_bias_model)) n_pps += 1;

  // pseudo-parameter for the Kn/Ks of each select-matrix category when
  // each smc has more than 1 parameter:
  const unsigned nsmc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_MATRIX));
  for(unsigned i(0);i<nsmc;++i){
    if(paramset_aa_select_matrix_cat_size(i)>1){
      n_pps += 1;
    }
  }

  return base_t::pseudo_param_size()+n_pps;
}



void
rate_gtor_nscodon_base::
set_pseudo_param(const unsigned param_no,
                 const smlfloat val) {

  if(param_no>=pseudo_param_size()){
    log_os << "ERROR:: Bad pseudoparam number\n";
    abort();
  } else if (param_no < base_t::pseudo_param_size()){
    base_t::set_pseudo_param(param_no,val);
  } else {
    // the pseudo-parameter is valid and handled in this class, so:
    const unsigned cat_no(param_no-base_t::pseudo_param_size());

    const bool is_cbias_pparam(is_cbias_pseudo_param(_data.copt.codon_bias_model));
    if(is_cbias_pparam && cat_no==0){
      // power-scale codon bias parameters
      const unsigned s(paramset_size(PARAMSET_CODON_BIAS));
      const unsigned st(paramset_start(PARAMSET_CODON_BIAS));
      for(unsigned i(st);i<(st+s);++i) {
        if(_is_train_param[i]) _param[i] = std::pow(_param[i],val);
      }
    } else {
      // scale dN/dS directly for each select matrix cat that has gt 1
      // param:
      unsigned dnds_scale_cat_no(cat_no);
      if(is_cbias_pparam) dnds_scale_cat_no--;
      unsigned n_pps(0);
      const unsigned nsmc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_MATRIX));
      for(unsigned i(0);i<nsmc;++i){
        if(paramset_aa_select_matrix_cat_size(i)>1){
          if(n_pps==dnds_scale_cat_no) {
            dnds_scale_cat_no=i;
            break;
          }
          n_pps += 1;
        }
      }

      assert(dnds_scale_cat_no<nsmc);
      {
        const smlfloat safe_val(std::abs(val));
        const unsigned s(paramset_aa_select_matrix_cat_size(dnds_scale_cat_no));
        const unsigned st(paramset_aa_select_matrix_cat_start(dnds_scale_cat_no));
        for(unsigned i(st);i<(st+s);++i) if(_is_train_param[i]) _param[i]*=safe_val;
      }
    }
  }
}




void
rate_gtor_nscodon_base::
nsaa_selection_matrix(irv_t<smlfloat> nsaa_selection[NSAA::SIZE*NSAA::SIZE],
                      const unsigned cat_no,
                      const unsigned branch_cat_set) const {

  for(unsigned i(0);i<NSAA::SIZE;++i){
    for(unsigned j(0);j<NSAA::SIZE;++j){
      if (i==j) {
        nsaa_selection[j+i*NSAA::SIZE] = 1.;
      } else {
        nsaa_selection[j+i*NSAA::SIZE] = cat_nsaa_param(static_cast<NSAA::index_t>(i),
                                                        static_cast<NSAA::index_t>(j),
                                                        cat_no,
                                                        branch_cat_set);
      }
    }
  }
}



void
rate_gtor_nscodon_base::
fix_input_param_internal(){

  base_t::fix_input_param_internal();
  {
    const CAT_PARAM_TYPE::index_t pt(CAT_PARAM_TYPE::SEL_STRENGTH);
    CAT_MIX_TYPE::index_t mt(CAT_MIX_TYPE::SITE);
    fix_input_param_cat_branch_paramset(PARAMSET_SITE_SELECT_CATS,pt,mt,
                                        cm().sel_strength_norm_set_size(mt),
                                        cm().sel_strength_norm_set_from_sel_strength_cat(mt));
    mt=CAT_MIX_TYPE::GROUP;
    fix_input_param_cat_branch_paramset(PARAMSET_GROUP_SELECT_CATS,pt,mt,
                                        cm().sel_strength_norm_set_size(mt),
                                        cm().sel_strength_norm_set_from_sel_strength_cat(mt));
  }

  if(_data.copt.codon_bias_model == CODON_BIAS_MODEL::SYNON_RATIO ||
     _data.copt.codon_bias_model == CODON_BIAS_MODEL::SYNON_NORM_RATIO ||
     _data.copt.codon_bias_model == CODON_BIAS_MODEL::SYNON_NORM_ABS){
    const unsigned pst(paramset_start(PARAMSET_CODON_BIAS));
    for(unsigned i(0);i<NSAA::SIZE;++i){
      const unsigned st(pst+_aacodon.aa_codon_start_offset[i]);
      const unsigned s(_aacodon.aa_codon_size[i]);
      if(s) {
        pdistro_norm_free_param(_param.begin()+st,
                                _param.begin()+st+s,
                                _is_train_param.begin()+st);
      }
    }
  }
}




unsigned
rate_gtor_nscodon_base::
paramset_size(const int p) const {
  if(p<base_t::PARAMSET_END){
    return base_t::paramset_size(p);
  }

  switch(static_cast<const paramset_t>(p)){
  case PARAMSET_SITE_SELECT_CATS :
    return cm().typed_cat_size(CAT_PARAM_TYPE::SEL_STRENGTH,CAT_MIX_TYPE::SITE);
  case PARAMSET_GROUP_SELECT_CATS :
    return cm().typed_cat_size(CAT_PARAM_TYPE::SEL_STRENGTH,CAT_MIX_TYPE::GROUP);
  case PARAMSET_AA :
    {
      unsigned psize(0);
      const unsigned nsmc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_MATRIX));
      for(unsigned i(0);i<nsmc;++i){
        psize += paramset_aa_select_matrix_cat_size(i);
      }
      return psize;
    }
  case PARAMSET_CODON_BIAS :
    {
      using namespace CODON_BIAS_MODEL;
      const index_t& cbm(_data.copt.codon_bias_model);
      if(cbm == SYNON_RATIO ||
         cbm == SYNON_NORM_RATIO ||
         cbm == SYNON_NORM_ABS) return NSCODON::SIZE;
    }
  default : return 0;
  }
}


// produce a human readable name for the aa selection parameters
//
void
rate_gtor_nscodon_base::
get_aa_param_index_name(const unsigned matrix_cat_no,
                        const unsigned matrix_param_no,
                        std::ostream& os) const {

  if(matrix_param_no==0){
    os << "aa_selection_matrix_" << matrix_cat_no << "_scale_param";
    return;
  }

  const unsigned sm_param_no(matrix_param_no-1);
  const SELECT_MODEL::index_t sm(select_model(matrix_cat_no));
  assert(sm_param_no<paramsize_select_model(sm));

  os << "aa_selection_matrix_" << matrix_cat_no;

  if       (sm == SELECT_MODEL::ASYMM){
    const unsigned aa_index(_aaex.asymm_paramtoaa_map[sm_param_no]);
    const unsigned from(aa_index/NSAA::SIZE);
    const unsigned to(aa_index%NSAA::SIZE);
    os << "_" << NSAA::syms[from] << "->" << NSAA::syms[to];
  } else if(sm == SELECT_MODEL::SYMM){
    const unsigned aa_index(_aaex.symm_paramtoaa_map[sm_param_no]);
    unsigned from(aa_index/NSAA::SIZE);
    unsigned to(aa_index%NSAA::SIZE);
    if(from>to) std::swap(from,to);
    os << "_" << NSAA::syms[from] << "<->" << NSAA::syms[to];
  } else if(sm == SELECT_MODEL::HP){
    const unsigned from(sm_param_no/NSAA_MAPS::HP::SIZE);
    const unsigned to(sm_param_no%NSAA_MAPS::HP::SIZE);
    os << "-" << NSAA_MAPS::HP::syms[from] << "->" << NSAA_MAPS::HP::syms[to];
  }

  os << "_param_" << sm_param_no;
}



void
rate_gtor_nscodon_base::
paramset_label(const paramset_t p,
               const unsigned n,
               std::ostream& os) const {

  switch(p){
  case PARAMSET_SITE_SELECT_CATS :
    os << "site_select_cats_" << n;
    break;
  case PARAMSET_GROUP_SELECT_CATS :
    os << "group_select_cats_" << n;
    break;
    //  case PARAMSET_GROUP_SELECT_CATS_PROB :
    //    os << "group_select_cats_prior_" << n;
    break;
  case PARAMSET_AA :
    {
      bool is_valid(false);
      unsigned matrix_cat_no(0);
      unsigned matrix_param_no(0);
      {
        const unsigned nsmc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_MATRIX));
        int ncopy(n);
        for(unsigned i(0);i<nsmc;++i){
          const int s(paramset_aa_select_matrix_cat_size(i));
          if(ncopy<s){
            if(ncopy>=0){
              matrix_cat_no=i;
              matrix_param_no=ncopy;
              is_valid=true;
            }
            break;
          }
          ncopy -= s;
        }
      }
      if(! is_valid) die("invalid paramset_label call");
      get_aa_param_index_name(matrix_cat_no,matrix_param_no,os);
    }
    break;
  case PARAMSET_CODON_BIAS :
    {
      const NSCODON::index_t c(_aacodon.nscodon_aa2normal_order_map[n]);
      os << "codon_bias_" << NSCODON::print(c,true);
    }
    break;
  default :
    os << "unknown_"  << n;
  }
}




const char* const section_id = "rate_gtor_nscodon_base";
const char* const end_label = "END";

const iliner il(section_id);



void
rate_gtor_nscodon_base::
store_state(std::ostream& os) const {

  base_t::store_state(os);

  os << section_id << " codon_bias_model " << _data.copt.codon_bias_model << "\n";

  const unsigned nsmc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_MATRIX));
  for(unsigned i(0);i<nsmc;++i){
    os << section_id << " select_model_" << i << " " << select_model(i) << "\n";
  }

  // write all paramsets:
  //
  for(unsigned p(base_t::PARAMSET_END);p<PARAMSET_END;++p){
    const unsigned st(paramset_start(static_cast<paramset_t>(p)));
    const unsigned s(paramset_size(static_cast<paramset_t>(p)));
    for(unsigned i(0);i<s;++i){
      os << section_id << " ";
      paramset_label(static_cast<paramset_t>(p),i,os);
      os << " " << _is_train_param[st+i] << " " << _param[st+i] << "\n";
    }
  }

  os << section_id << " " << end_label << "\n";
}



void
rate_gtor_nscodon_base::
load_state_internal(std::istream& is) {

  base_t::load_state_internal(is);

  il.advance(is,"codon_bias_model");
  _data.copt.codon_bias_model=static_cast<CODON_BIAS_MODEL::index_t>(read_int(is));

  const unsigned nsmc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_MATRIX));
  for(unsigned i(0);i<nsmc;++i){
    std::ostringstream oss;
    oss << "select_model_" << i;
    il.advance(is,oss.str().c_str());
    const SELECT_MODEL::index_t sm(static_cast<SELECT_MODEL::index_t>(read_int(is)));
    _data.copt.set_select_model(i,sm);
  }

  init_class_storage();

  // read all paramsets:
  //
  for(unsigned p(base_t::PARAMSET_END);p<PARAMSET_END;++p){
    const unsigned st=paramset_start(static_cast<paramset_t>(p));
    const unsigned s=paramset_size(static_cast<paramset_t>(p));
    for(unsigned i(0);i<s;++i){
      il.advance(is);
      bool tmpb;
      is >> tmpb >> _param[st+i];
      _is_train_param[st+i] = tmpb;
    }
  }

  il.advance(is,end_label);
}


NSAA_MAPS::converter aacon;

// get aa parameter index from aa codes:
// return -1 for fixed values
std::pair<unsigned,AA_PARAM::index_t>
rate_gtor_nscodon_base::
nsaa_param_index(const NSAA::index_t from,
                 const NSAA::index_t to,
                 const unsigned select_matrix_cat) const {

  assert( to != from );
  const unsigned aa_index(to+from*NSAA::SIZE);

  if(! _aaex.is_1nuc_aa_exchange[aa_index]) return std::make_pair(0,AA_PARAM::zero);;

  using namespace SELECT_MODEL;
  const index_t s(select_model(select_matrix_cat));
  unsigned aa_offset;
  if       (s == NONE || s == FROM || s == TO || s == FROM_TO){
    return std::make_pair(0,AA_PARAM::one);
  } else if(s == SINGLE){
    aa_offset = 0;
  } else if(s == HP){
    const NSAA_MAPS::HP::index_t hpf(aacon.aa_to_hp(from));
    const NSAA_MAPS::HP::index_t hpt(aacon.aa_to_hp(to));
    aa_offset = hpf*NSAA_MAPS::HP::SIZE+hpt;
  } else if(s == SYMM || s == FROM_SYMM || s == SYMM_TO || s == FROM_SYMM_TO){
    aa_offset = _aaex.symm_aatoparam_map[aa_index];
  } else if(s == ASYMM){
    aa_offset = _aaex.asymm_aatoparam_map[aa_index];
  } else { die("Invalid Selection model"); }

  const unsigned cat_param_start(paramset_aa_select_matrix_cat_start(select_matrix_cat));

  unsigned param_index(cat_param_start+aa_offset);

  if( is_sm_scaled(s)) param_index += 1; // leave room for the scale param

  assert(param_index < paramset_start(PARAMSET_AA)+paramset_size(PARAMSET_AA));
  return std::make_pair(param_index,AA_PARAM::normal);
}



// get aa parameter from aa codes:
irv_t<smlfloat>
rate_gtor_nscodon_base::
nsaa_param(const NSAA::index_t from,
           const NSAA::index_t to,
           const unsigned select_matrix_cat) const {


  irv_t<smlfloat> val;
  {
    std::pair<unsigned,AA_PARAM::index_t> p = nsaa_param_index(from,to,select_matrix_cat);
    if( p.second == AA_PARAM::one ) {
      val=irv_t<smlfloat>(1.);
    } else if ( p.second == AA_PARAM::zero ) {
      val=irv_t<smlfloat>(0.);
    } else {
      val=irv_t<smlfloat>(_param[p.first],param_variance(p.first));
    }
  }

  {
    using namespace SELECT_MODEL;
    index_t s(select_model(select_matrix_cat));
    const unsigned st(paramset_aa_select_matrix_cat_start(select_matrix_cat));

    // factor in other parameters:
    if       (s == FROM){
      const unsigned aafactor_index(st+1+static_cast<unsigned>(from));
      val = val * irv_t<smlfloat>(_param[aafactor_index],param_variance(aafactor_index));
    } else if(s == TO){
      const unsigned aafactor_index(st+1+static_cast<unsigned>(to));
      val = val * irv_t<smlfloat>(_param[aafactor_index],param_variance(aafactor_index));
    } else if(s == FROM_TO){
      const unsigned aafactor_index1(st+1+static_cast<unsigned>(from));
      const unsigned aafactor_index2(st+1+NSAA::SIZE+static_cast<unsigned>(to));
      val = val * irv_t<smlfloat>(_param[aafactor_index1],param_variance(aafactor_index1));
      val = val * irv_t<smlfloat>(_param[aafactor_index2],param_variance(aafactor_index2));
    } else if(s == FROM_SYMM){
      const unsigned aafactor_index(st+_aaex.symm_size+1+static_cast<unsigned>(from));
      val = val * irv_t<smlfloat>(_param[aafactor_index],param_variance(aafactor_index));
    } else if(s == SYMM_TO){
      const unsigned aafactor_index(st+_aaex.symm_size+1+static_cast<unsigned>(to));
      val = val * irv_t<smlfloat>(_param[aafactor_index],param_variance(aafactor_index));
    } else if(s == FROM_SYMM_TO){
      const unsigned aafactor_index1(st+_aaex.symm_size+1+static_cast<unsigned>(from));
      const unsigned aafactor_index2(st+_aaex.symm_size+1+NSAA::SIZE+static_cast<unsigned>(to));
      val = val * irv_t<smlfloat>(_param[aafactor_index1],param_variance(aafactor_index1));
      val = val * irv_t<smlfloat>(_param[aafactor_index2],param_variance(aafactor_index2));
    }

    // scale parameter
    if( is_sm_scaled(s) ){
      val = val * irv_t<smlfloat>(_param[st],param_variance(st));
    }

    return val;
  }
}



irv_t<smlfloat>
rate_gtor_nscodon_base::
cat_nsaa_param(NSAA::index_t from,
               NSAA::index_t to,
               const unsigned cat_no,
               const unsigned branch_cat_set
               //const SCALE_TYPE::index_t st
               ) const {

  if(from==to) return 1.;

  const unsigned scm(cm().typed_cat_no_from_cat_no_y_branch_cat_set(cat_no,branch_cat_set,CAT_PARAM_TYPE::SEL_MATRIX));
  irv_t<smlfloat> val(nsaa_param(from,to,scm));

  {
    const irv_t<smlfloat> sfactor = get_paramset_member(PARAMSET_SITE_SELECT_CATS,
                                                        cm().typed_cat_no_from_cat_no_y_branch_cat_set(cat_no,
                                                                                                       branch_cat_set,
                                                                                                       CAT_PARAM_TYPE::SEL_STRENGTH,
                                                                                                       CAT_MIX_TYPE::SITE));
    // \todo scale_type
    val *= sfactor;
  }
  {
    const irv_t<smlfloat> sfactor = get_paramset_member(PARAMSET_GROUP_SELECT_CATS,
                                                        cm().typed_cat_no_from_cat_no_y_branch_cat_set(cat_no,
                                                                                                       branch_cat_set,
                                                                                                       CAT_PARAM_TYPE::SEL_STRENGTH,
                                                                                                       CAT_MIX_TYPE::GROUP));
    val *= sfactor;
  }

  return val;
}



smlfloat
rate_gtor_nscodon_base::
codon_bias_factor(const NSCODON::index_t codon) const {
  const unsigned st(paramset_start(PARAMSET_CODON_BIAS));

  using namespace CODON_BIAS_MODEL;
  const index_t& cbm(_data.copt.codon_bias_model);

  if (cbm==SYNON_RATIO || cbm==SYNON_NORM_RATIO || cbm==SYNON_NORM_ABS){

    // get codons reindexed into continuous aa blocks
    const unsigned reindexed_codon(_aacodon.nscodon_normal2aa_order_map[codon]);
    const NSAA::index_t aa(codon_trans_known(codon));

    if(cbm == SYNON_NORM_ABS){
      const unsigned st2(_aacodon.aa_codon_start_offset[aa]);
      const unsigned s2(_aacodon.aa_codon_size[aa]);
      const std::vector<smlfloat>::const_iterator max_value_in_aa(std::max_element(_param.begin()+st+st2,_param.begin()+st+st2+s2));
      return _param[st+reindexed_codon]/(*max_value_in_aa);
    } else {
      // the codon_size multiplier here is the completion of the
      // normalization step, to make the expectation of each aa's codons
      // equal to 1.-> obviously it would be better to complete
      // normalization in the normalization function
      return _param[st+reindexed_codon]*static_cast<smlfloat>(_aacodon.aa_codon_size[aa]);
    }
  } else {
    die("invalid codon_bias_factor call");
  }
}




smlfloat
rate_gtor_nscodon_base::
codon_bias(const NSCODON::index_t from,
           const NSCODON::index_t to) const {

  using namespace CODON_BIAS_MODEL;

  const index_t& cbm(_data.copt.codon_bias_model);

  if       (cbm == SYNON_NORM_RATIO){
    return codon_bias_factor(to)/codon_bias_factor(from);
  } else if(cbm == SYNON_RATIO){
    if(codon_trans_known(to)!=codon_trans_known(from)) return 1.;
    return codon_bias_factor(to)/codon_bias_factor(from);
  } else if(cbm == NONE){
    return 1.;
  } else {
    die("unknown codon bias model type");
  }
}



smlfloat
rate_gtor_nscodon_base::
param_penalty() const {
  static const smlfloat paa_penalty_val(50.);
  static const smlfloat paa_penalty_factor(0.0001);
  const unsigned st(paramset_start(PARAMSET_AA));
  const unsigned s(paramset_size(PARAMSET_AA));
  smlfloat pen(0.);
  for(unsigned i(st);i<(st+s);++i) {
    if(_param[i]>paa_penalty_val) {
      const smlfloat x(1.+(_param[i]-paa_penalty_val)*paa_penalty_factor);
      pen += (x*x)-1.;
    }
  }
  return pen;
}



void
rate_gtor_nscodon_base::
state_pdistro_to_nuc_pdistro(const prob_t* state_pdistro,
                             prob_t* nuc_pdistro) const {
  prob_t nscodon_pdistro[NSCODON::SIZE];
  state_pdistro_to_nscodon_pdistro(state_pdistro,nscodon_pdistro);
  nscodon_pdf_2_nuc_pdf(nscodon_pdistro,nuc_pdistro);
}


void
rate_gtor_nscodon_base::
state_pdistro_to_nsaa_pdistro(const prob_t* state_pdistro,
                              prob_t* nsaa_pdistro) const {
  prob_t nscodon_pdistro[NSCODON::SIZE];
  state_pdistro_to_nscodon_pdistro(state_pdistro,nscodon_pdistro);
  nscodon_pdf_2_nsaa_pdf(nscodon_pdistro,nsaa_pdistro);
}


void
rate_gtor_nscodon_base::
bg_pdistro_nuc_pos_update(){
  const unsigned bcn(bg_cat_size());

  _data.bg_pdistro_nuc_pos_bg_cat.init(bcn,CODON::BASE_SIZE,NUC::SIZE);
  for(unsigned c(0);c<bcn;++c){
    for(unsigned pos(0);pos<CODON::BASE_SIZE;++pos){
      nscodon_pdf_2_nuc_pdf_pos(bg_pdistro_nscodon_bg_cat(c),pos,bg_pdistro_nuc_pos_bg_cat(c,pos));
    }
  }
}



static
void
bg_pdistro_nuc_pos_conditioned_update(const SITE_MODEL::index_t sm,
                                      const prob_t* bg_pdistro,
                                      prob_t** bg_pdistro_nuc_pos_conditioned_on_3p,
                                      prob_t** bg_pdistro_nuc_pos_conditioned_on_5p){

  const unsigned base_size(SITE_MODEL::base_size(sm));
  unsigned codon_pos[SITE_MODEL::MAX_BASE_SIZE];
  SITE_MODEL::codon_position(sm,codon_pos);

  for(unsigned pos(0);pos<CODON::BASE_SIZE;++pos){
    unsigned base_no(0);
    {
      bool is_base_no_3p(false);
      unsigned base_no_3p(0);

      for(unsigned b(0);b<base_size;++b){
        if(codon_pos[b] == pos){
          base_no=b;
          if((b+1)<base_size){
            base_no_3p=b+1;
            is_base_no_3p=true;
            break;
          }
        }
      }

      get_dependent_nuc_pos_distro_from_site_model_distro(sm,bg_pdistro,base_no,is_base_no_3p,base_no_3p,
                                                          bg_pdistro_nuc_pos_conditioned_on_3p[pos]);
    }
    {
      bool is_base_no_5p(false);
      unsigned base_no_5p(0);

      for(unsigned b(0);b<base_size;++b){
        if(codon_pos[b] == pos){
          base_no=b;
          if(b>0){
            base_no_5p=b-1;
            is_base_no_5p=true;
            break;
          }
        }
      }

      get_dependent_nuc_pos_distro_from_site_model_distro(sm,bg_pdistro,base_no,is_base_no_5p,base_no_5p,
                                                          bg_pdistro_nuc_pos_conditioned_on_5p[pos]);
    }
  }
}




void
rate_gtor_nscodon_base::
bg_pdistro_nuc_pos_conditioned_update(){
  const SITE_MODEL::index_t sm(site_model());
  const unsigned bcn(bg_cat_size());
  _data.bg_pdistro_nuc_pos_conditioned_on_3p_bg_cat.init(bcn,CODON::BASE_SIZE,NUC::SIZE*NUC::SIZE);
  _data.bg_pdistro_nuc_pos_conditioned_on_5p_bg_cat.init(bcn,CODON::BASE_SIZE,NUC::SIZE*NUC::SIZE);
  for(unsigned c(0);c<bcn;++c){
    ::bg_pdistro_nuc_pos_conditioned_update(sm,bg_pdistro_bg_cat(c),
                                            _data.bg_pdistro_nuc_pos_conditioned_on_3p_bg_cat[c],
                                            _data.bg_pdistro_nuc_pos_conditioned_on_5p_bg_cat[c]);
  }


  //#define POS_DEP_DEBUG
#ifdef POS_DEP_DEBUG
  log_os << "codon pos nuc distros:\n";
  for(unsigned pos(0);pos<CODON::BASE_SIZE;++pos){
    log_os << "pos: " << pos << " :";
    for(unsigned i(0);i<NUC::SIZE;++i){
      log_os << " " << bg_pdistro_nuc_pos_cat(0,pos)[i];
    }
    log_os << "\n";
  }
  log_os << "\n";

  log_os << std::setprecision(8);
  log_os << "codon pos dep_5p nuc distros:\n";
  for(unsigned pos(0);pos<CODON::BASE_SIZE;++pos){
    for(unsigned i(0);i<NUC::SIZE;++i){
      log_os << "pos: " << pos << " [" << NUC::syms[i] <<"]:";
      for(unsigned j(0);j<NUC::SIZE;++j){
        log_os << " " << bg_pdistro_nuc_pos_conditioned_on_5p_cat(0,pos,static_cast<NUC::index_t>(i))[j];
      }
      log_os << "\n";
    }
  }
  log_os << "\n";

  log_os << "codon pos dep_3p nuc distros:\n";
  for(unsigned pos(0);pos<CODON::BASE_SIZE;++pos){
    for(unsigned i(0);i<NUC::SIZE;++i){
      log_os << "pos: " << pos << " [" << NUC::syms[i] <<"]:";
      for(unsigned j(0);j<NUC::SIZE;++j){
        log_os << " " << bg_pdistro_nuc_pos_conditioned_on_3p_cat(0,pos,static_cast<NUC::index_t>(i))[j];
      }
      log_os << "\n";
    }
  }
#endif

}


void
rate_gtor_nscodon_base::
stat_pdistro_nscodon(prob_t* nscodon_sd,
                     const rates_func_options& opt) const {
  prob_t* sd(new prob_t[state_size()]);
  get_stationary_pdistro(sd,*this,opt);
  state_pdistro_to_nscodon_pdistro(sd,nscodon_sd);
  delete sd;
}



void
rate_gtor_nscodon_base::
stat_pdistro_nsaa(prob_t* nsaa_sd,
                  const rates_func_options& opt) const {
  prob_t* sd(new prob_t[state_size()]);
  get_stationary_pdistro(sd,*this,opt);
  state_pdistro_to_nsaa_pdistro(sd,nsaa_sd);
  delete sd;
}



void
rate_gtor_nscodon_base::
reset_param_internal(const PARAM_INIT_TYPE::index_t pinit){

  base_t::reset_param_internal(pinit);

  if(pinit == PARAM_INIT_TYPE::RANDOM){
    static const smlfloat select_cats_min(0.1);
    static const smlfloat select_cats_range(2.);

    {
      const unsigned st = paramset_start(PARAMSET_SITE_SELECT_CATS);
      const unsigned s = paramset_size(PARAMSET_SITE_SELECT_CATS);
      for(unsigned i(st);i<(st+s);++i){
        if(_is_train_param[i]) { _param[i] = random_uniform()*select_cats_range+select_cats_min; }
        else                   { _param[i] = 1.; }
      }
    }

    {
      const unsigned st = paramset_start(PARAMSET_GROUP_SELECT_CATS);
      const unsigned s = paramset_size(PARAMSET_GROUP_SELECT_CATS);
      for(unsigned i(st);i<(st+s);++i){
        if(_is_train_param[i]) { _param[i] = random_uniform()*select_cats_range+select_cats_min; }
        else                   { _param[i] = 1.; }
      }
    }

    {
      /// \todo make these values reasonable for newer select models:
      static const smlfloat prate_min(0.05);
      static const smlfloat prate_max(1.5);
      static const smlfloat praterange(prate_max-prate_min);

      const unsigned staa = paramset_start(PARAMSET_AA);
      const unsigned saa = paramset_size(PARAMSET_AA);
      for(unsigned i(staa);i<(staa+saa);++i){ _param[i] = random_uniform()*praterange+prate_min; }
    }

    {
      static const smlfloat codon_bias_min(0.1);
      static const smlfloat codon_bias_range(0.9);

      const unsigned st(paramset_start(PARAMSET_CODON_BIAS));
      const unsigned s(paramset_size(PARAMSET_CODON_BIAS));
      for(unsigned i(st);i<(st+s);++i){
        if(_is_train_param[i]) { _param[i] = random_uniform()*codon_bias_range+codon_bias_min; }
        else                   { _param[i] = 1.; }
      }
    }

  } else if(pinit == PARAM_INIT_TYPE::START){
    // don't setup initial category parameters on a slope if they are
    // tied to assigned data
    if(cm().assigned_data_set_size() > 1) die("can't init params with multi adsets");
    set_paramset_slope(PARAMSET_SITE_SELECT_CATS);
    set_paramset_slope(PARAMSET_GROUP_SELECT_CATS);

    {
      using namespace SELECT_MODEL;

      static const smlfloat start_select_base(0.15);
      smlfloat start_select(start_select_base);

      const unsigned nsmc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_MATRIX));
      for(unsigned i(0);i<nsmc;++i){
        const index_t sm(select_model(i));
        if(sm==FROM_SYMM_TO){
          static const smlfloat onethird(1./3.);
          start_select=std::pow(start_select_base,onethird);
        } else if(sm==FROM_TO || sm== FROM_SYMM || sm==SYMM_TO){
          start_select=std::sqrt(start_select_base);
        }

        const unsigned st(paramset_aa_select_matrix_cat_start(i));
        const unsigned si(paramset_aa_select_matrix_cat_size(i));
        for(unsigned j(st);j<(st+si);++j){ _param[j] = start_select; }
      }
    }
    set_paramset_val(PARAMSET_CODON_BIAS,1.);
  }

  const unsigned nsmc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_MATRIX));
  for(unsigned i(0);i<nsmc;++i){
    using namespace SELECT_MODEL;

    const index_t sm(select_model(i));
    const unsigned st(paramset_aa_select_matrix_cat_start(i));
    if(!(sm==NONE||sm==SINGLE)){ _param[st] = 1.;}
  }
}


// i686 gcc 4.1 has lots of (harmless) link warnings when this is left inline:
//
void
rate_gtor_nscodon_base::
bg_pdistro_update_internal(){
  base_t::bg_pdistro_update_internal();
  bg_pdistro_nscodon_update();
  bg_pdistro_nuc_pos_update();
  bg_pdistro_nuc_pos_conditioned_update();
  _data.is_bg_valid=true;
}



smlfloat
rate_gtor_nscodon_base::
cat_nuc_context_mut_rate(const NUC::index_t from0,
                         const NUC::index_t from1,
                         const NUC::index_t from2,
                         const NUC::index_t to,
                         const unsigned cat_no,
                         const unsigned branch_cat_set,
                         const prob_t* nuc_bg0,
                         const prob_t* nuc_bg2,
                         const NSCODON::index_t to_state) const {

  smlfloat val(base_t::cat_nuc_context_mut_rate(from0,from1,from2,to,cat_no,branch_cat_set,nuc_bg0,nuc_bg2,true));

  const unsigned mmc(cm().typed_cat_no_from_cat_no_y_branch_cat_set(cat_no,branch_cat_set,CAT_PARAM_TYPE::MUT_MODEL));
  if(rate_model_nuc(mmc)==RATE_MODEL_NUC::GY94){
    if(model_type()!=RATE_GTOR_MODEL::CODON){
      die("GY94 rates are valid for codon model only");
    }
    // scale by nscodon size to keep tree branches in valid range.
    val *= bg_pdistro_nscodon_cat(cat_no)[to_state]*static_cast<smlfloat>(NSCODON::SIZE);
  }
  return val;
}

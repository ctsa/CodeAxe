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

// $Id: rate_gtor_nuc_base.cc 1213 2008-05-21 23:37:51Z ctsa $

/// \file

#include "cat_manager.h"
#include "context_group.h"
#include "rate_gtor_nuc_base.h"
#include "rate_gtor_util.h"
#include "rate_gtor_nuc_base_util.h"
#include "simple_util.h"
#include "substk_exception.h"
#include "util/bio/bioseq_util_flux_rate_conversions.h"
#include "util/general/die.h"
#include "util/general/io_util.h"
#include "util/general/log.h"
#include "util/math/indy_random_var_io.h"
#include "util/math/matrix_util_io.h"
#include "util/math/polynom_util.h"
#include "util/math/prob_util_io.h"
#include "util/math/random_util.h"

#include <cassert>

#include <ostream>
#include <utility>
#include <vector>


static
unsigned
get_rate_model_param_size(const RATE_MODEL_NUC::index_t i){

  using namespace RATE_MODEL_NUC;

  switch(i){
  case JC69:
  case F81:    return 1;
  case K80:
  case HKY85:
  case GY94:   return 2;
  case REV:    return (off_diag_size(NUC::SIZE)/2);
  case NONREV: return off_diag_size(NUC::SIZE);
  default: die("Unsupported nuc rate model");
  }
}



static
unsigned
context_adjust_nuc_param_size(const CONTEXT_MODEL_NUC::index_t i,
                              const unsigned s){

  using namespace CONTEXT_MODEL_NUC;

  switch(i){
  case INDY:             return s;
  case CPG_1TI:          return s+1;
  case CPG_ONLY:
  case TPA_ONLY:
  case CPG_2TI:          return s+2;
  case CPG_NONREV:       return s+2*(NUC::SIZE-1);
  case PRE_DOUBLET:
  case POST_DOUBLET:     return s*NUC::SIZE;
  case FACTORED_TRIPLET: return s*2*NUC::SIZE;
  case TRIPLET:          return s*NUC::SIZE*NUC::SIZE;
  default: die("Unsupported nuc context model");
  }
}



unsigned
rate_gtor_nuc_base::
get_rate_context_model_param_size(const RATE_MODEL_NUC::index_t rm,
                                  const CONTEXT_MODEL_NUC::index_t cm) {
  unsigned model_param_size(get_rate_model_param_size(rm));
  return context_adjust_nuc_param_size(cm,model_param_size);
}



void
rate_gtor_nuc_base::
init_class_storage(){

  const unsigned nmmc(cm().typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL));

  _data.nopt.set_n_site_mutation_model_cats(nmmc);

  {
    bool is_all_indy(true);
    for(unsigned i(0);i<nmmc;++i){
      if(context_model_nuc(i)!=CONTEXT_MODEL_NUC::INDY){
        is_all_indy=false;
      }
    }

    if(is_all_indy) _data.nopt.is_use_edge_correction = false;
  }

  resize(param_end());

  {
    const unsigned st(base_t::param_end());
    const unsigned sts(param_end());
    for(unsigned i(st);i<sts;++i){
      _is_free_block_start[i] = false;
      _is_free_block_stop[i] = false;
    }
  }

  /// \todo update param_object and this code to allow non-continuous normalization blocks
  ///       also, I don't think this dependent block is even accurate anymore???
  set_paramset_dependent_block(PARAMSET_SITE_RATE_CATS);
  set_paramset_dependent_block(PARAMSET_GROUP_RATE_CATS);

  // independent blocks of mutation model parameters can be
  // discontinuous by the catmodel grammer, but this case isn't
  // supported yet by the param_object indy block transformer
  //
  for(unsigned i(0);i<nmmc;++i){
    const unsigned st(paramset_mutation_model_cat_start(i));
    const unsigned s(paramset_mutation_model_cat_size(i));

    const unsigned set_no(cm().mut_model_norm_set_from_mut_model_cat()[i]);

    const bool is_first(i==0);
    unsigned last_set_no(set_no);
    if(! is_first){
      last_set_no=cm().mut_model_norm_set_from_mut_model_cat()[i-1];
    }
    if(is_first || last_set_no != set_no){
      if( (! is_first) && last_set_no+1 != set_no){
        pass_away("can't currently handle discontinuous normalization blocks");
      }
      _is_free_block_start[st] = true;
    }

    const bool is_last((i+1)==nmmc);
    unsigned next_set_no(set_no);
    if(! is_last){
      next_set_no=cm().mut_model_norm_set_from_mut_model_cat()[i+1];
    }
    if(is_last || next_set_no != set_no){
      if( (! is_last) && next_set_no != set_no+1){
        pass_away("can't currently handle discontinuous normalization blocks");
      }
      _is_free_block_stop[st+s-1] = true;
    }
  }
}



void
rate_gtor_nuc_base::
lock_1param_norm_sets(const int paramset,
                      const CAT_PARAM_TYPE::index_t pt,
                      const CAT_MIX_TYPE::index_t mt,
                      const unsigned n_sets,
                      const std::vector<unsigned>::const_iterator set_from_cat){

  const unsigned st(paramset_start(paramset));
  const unsigned t_cat_size(cm().typed_cat_size(pt,mt));

  for(unsigned j(0);j<n_sets;++j){
    unsigned block_count(0);
    for(unsigned i(0);i<t_cat_size;++i){
      if(set_from_cat[i]==j) block_count++;
    }
    if(block_count==1){
      for(unsigned i(0);i<t_cat_size;++i){
        if(set_from_cat[i]==j) _is_train_param[st+i]=false;
      }
    }
  }
}



void
rate_gtor_nuc_base::
init_param_locks(const rate_gtor_nuc_lock_options& lopt){

  // lock parameters which won't be trained
  {
    CAT_PARAM_TYPE::index_t pt(CAT_PARAM_TYPE::MUT_RATE);
    if(lopt.is_lock_site_rate_cats){
      set_paramset_train(PARAMSET_SITE_RATE_CATS,false);
    } else {
      CAT_MIX_TYPE::index_t mt(CAT_MIX_TYPE::SITE);
      lock_1param_norm_sets(PARAMSET_SITE_RATE_CATS,pt,mt,
                            cm().mut_rate_norm_set_size(mt),
                            cm().mut_rate_norm_set_from_mut_rate_cat(mt));
    }

    if(lopt.is_lock_group_rate_cats){
      set_paramset_train(PARAMSET_GROUP_RATE_CATS,false);
    } else {
      CAT_MIX_TYPE::index_t mt(CAT_MIX_TYPE::GROUP);
      lock_1param_norm_sets(PARAMSET_GROUP_RATE_CATS,pt,mt,
                            cm().mut_rate_norm_set_size(mt),
                            cm().mut_rate_norm_set_from_mut_rate_cat(mt));
    }
  }

  /// \todo figure out a good solution to factored triplet
  /// overdetermined parameter problem, as described below:
  //
  // A  different   type  of  parameter  dependence   exists  for  the
  // factored-triplet context  model. Other overdetermined  sets (such
  // as prob distros) can be  arbitrarily scaled, but we cannot put an
  // absolute  lock on  a single  parameter. For  the factored-triplet
  // model  the opposite is  true for  the parameters  associated with
  // each  from->to mutation:  we cannot  scale, but  we can  lock one
  // parameter (and still achieve the  full range of values). For this
  // reason, we arbitrarily lock the first free parameter for each set
  // here.  This also  means  that a  uh...   BUT HOW  WILL THIS  PLAY
  // TOGETHER  WITH  FREEBLOCK  scaling??????  I think  it  will  play
  // fine... if we consider only one set if the first member is locked
  // to  any constant,  the  group can  still  be scaled  to have  any
  // average rate.
  //
  // really sloppy solution -- produces valid ci's at least, but makes
  // minimization very slow, b/c of a) really bad start parameters, and
  // b) bad conj-dir minimization
  //

#if 0
  static const unsigned ft_block_size(NUC::SIZE+NUC::SIZE);
  const unsigned nmmc(cat_size_mutation_model());
  for(unsigned i(0);i<nmmc;++i){
    if(context_model_nuc(i)!=CONTEXT_MODEL_NUC::FACTORED_TRIPLET) continue;

    const unsigned st(paramset_mutation_model_cat_start(i));
    const unsigned s(paramset_mutation_model_cat_size(i));

    const unsigned n_blocks(s/ft_block_size);
    for(unsigned mut(0);mut<n_blocks;++mut){
      for(unsigned j(0);j<ft_block_size;++j){
        const unsigned param_no(st+mut*ft_block_size+j);
        if(_is_train_param[param_no]){
          _param[param_no]=1.;
          _is_train_param[param_no]=false;
          break;
        }
      }
    }
  }
#endif

#if 0
  // disabled this, b/c can't find a way to set val correctly first
  { // if jukes rate model, then lock nuc param
    const unsigned st(paramset_start(PARAMSET_NUC));
    const unsigned s(paramset_size(PARAMSET_NUC));
    if(s==1){
      _is_train_param[st]=false;
    }
  }
#endif
}




void
rate_gtor_nuc_base::
rates(smlfloat* rts,
      const rates_func_options& opt) const {
  rates_nuc_context(*this,rts,opt);
}



smlfloat
rate_gtor_nuc_base::
rate_scale_cat_no_y_branch_id(const unsigned cat_no,
                              const unsigned branch_id,
                              const trainable_mode::index_t tm) const {

  const prob_t* bgn(bg_pdistro_nuc_cat(cat_no));
  const unsigned branch_cat_set(cm().get_branch_cat_set(cat_no,branch_id));

  smlfloat sum(0.);
  for(unsigned n1(0);n1<NUC::SIZE;++n1){
    const NUC::index_t n1n(static_cast<NUC::index_t>(n1));
    for(unsigned n2(0);n2<NUC::SIZE;++n2){
      if(n1 == n2) continue;
      const NUC::index_t n2n(static_cast<NUC::index_t>(n2));
      const smlfloat val(nuc_context_mut_rate(NUC::N,n1n,
                                              NUC::N,n2n,
                                              cat_no,branch_cat_set,bgn,bgn,tm));

      sum += bgn[n1]*val;
    }
  }

  return sum;
}



void
rate_gtor_nuc_base::
stat_pdistro_nuc(prob_t* nuc_sd,
                 const rates_func_options& opt) const {
  prob_t* sd(new prob_t[state_size()]);
  get_stationary_pdistro(sd,*this,opt);
  state_pdistro_to_nuc_pdistro(sd,nuc_sd);
  delete [] sd;
}



void
rate_gtor_nuc_base::
fix_input_param_cat_branch_paramset(const int paramset,
                                    const CAT_PARAM_TYPE::index_t pt,
                                    const CAT_MIX_TYPE::index_t mt,
                                    const unsigned n_sets,
                                    const std::vector<unsigned>::const_iterator& set_from_cat){

  // run uncompensated normalization on free parameters in cat paramset
  const unsigned st(paramset_start(paramset));
  const unsigned s(paramset_size(paramset));

  assert(s==cm().typed_cat_size(pt,mt));

  const unsigned n_branches(cm().branch_size());
  const unsigned n_cats(cm().cat_size());
  simple_array<prob_t> cp(n_cats);
  cm().cat_pdistro(cp.ptr());

  const unsigned n_tcats(s);
  simple_array<prob_t> tcp(n_tcats,0.);

  for(unsigned c(0);c<n_cats;++c){
    const prob_t p(cp[c]);

    for(unsigned b(0);b<n_branches;++b){
      const unsigned tc(cm().typed_cat_no_from_cat_no_y_branch_id(c,b,pt,mt));
      tcp[tc] += p;
    }
  }

  for(unsigned set(0);set<n_sets;++set){
    simple_array<prob_t> tcp_set(tcp);
    simple_array<bool> is_train_tc_set(n_tcats);
    for(unsigned i(0);i<n_tcats;++i){
      is_train_tc_set[i] = _is_train_param[st+i];
    }

    for(unsigned tc(0);tc<n_tcats;++tc){
      const bool is_set(set==set_from_cat[tc]);
      if(!is_set) tcp_set[tc] = 0.;
      is_train_tc_set[tc] = is_train_tc_set[tc] && is_set;
    }

    static const smlfloat expect(1.);

    const smlfloat scale =
      array_scale_expectation_selected_only(_param.begin()+st,
                                            _param.begin()+st+s,
                                            tcp_set.ptr(),
                                            is_train_tc_set.ptr(),
                                            expect);
    if(scale<0) {
      throw substk_exception("rate_gtor_nuc_base.fix_input_param_cat_paramset(): locked parameter is out of range");
    }
  }
}



void
rate_gtor_nuc_base::
normalize_mut_rate_params(){
  const unsigned n_sets(cm().mut_model_norm_set_size());
  const unsigned n_cats(cm().cat_size());
  const unsigned n_branches(cm().branch_size());

  std::vector<smlfloat> set_rate_all(n_sets,0.);
  std::vector<smlfloat> set_rate_full(n_sets,0.);
  std::vector<smlfloat> set_rate_half(n_sets,0.);
  std::vector<prob_t> set_prob(n_sets,0.);

  simple_array<prob_t> cp(n_cats);
  cm().cat_pdistro(cp.ptr());

  // total normalization components for each mut model cat set:
  //
  for(unsigned c(0);c<n_cats;++c){
    const prob_t p(cp[c]/static_cast<prob_t>(n_branches));

    for(unsigned b(0);b<n_branches;++b){
      const unsigned mmc(cm().typed_cat_no_from_cat_no_y_branch_id(c,b,CAT_PARAM_TYPE::MUT_MODEL));
      const unsigned set(cm().mut_model_norm_set_from_mut_model_cat()[mmc]);

      // take into consideration which parameters are locked to
      // calculate normalization scale:
      //
      const smlfloat r_all(rate_scale_cat_no_y_branch_id(c,b,trainable_mode::NONE));
      const smlfloat r_full(rate_scale_cat_no_y_branch_id(c,b,trainable_mode::FULL));
      const smlfloat r_half(rate_scale_cat_no_y_branch_id(c,b,trainable_mode::HALF_JOINT));

      set_rate_all[set] += r_all*p;
      set_rate_full[set] += r_full*p;
      set_rate_half[set] += r_half*p;
      set_prob[set] += p;
    }
  }

  std::vector<smlfloat> mut_model_scale(n_sets,0.);

  static const smlfloat expect(1.);

  for(unsigned set(0);set<n_sets;++set){
    set_rate_all[set] /= set_prob[set];
    set_rate_full[set] /= set_prob[set];
    set_rate_half[set] /= set_prob[set];

    const smlfloat r_lock(set_rate_all[set] - ( set_rate_full[set] + set_rate_half[set] ));

    smlfloat s1,s2;
    solve_quadratic(set_rate_full[set],set_rate_half[set],r_lock-expect,s1,s2);

    mut_model_scale[set] = (s1*s1);
  }

  const unsigned nmmc(cm().typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL));
  for(unsigned mmc(0);mmc<nmmc;++mmc){
    const unsigned set(cm().mut_model_norm_set_from_mut_model_cat()[mmc]);

    scale_nuc_mut_rate_trainable_param(mmc,mut_model_scale[set]);
  }


  { /// finished scaling sets, so now run validation:

    std::vector<smlfloat> expect_scale(n_sets,0.);

    for(unsigned c(0);c<n_cats;++c){
      const prob_t p(cp[c]/static_cast<prob_t>(n_branches));

      for(unsigned b(0);b<n_branches;++b){
        const unsigned mmc(cm().typed_cat_no_from_cat_no_y_branch_id(c,b,CAT_PARAM_TYPE::MUT_MODEL));
        const unsigned set(cm().mut_model_norm_set_from_mut_model_cat()[mmc]);

        expect_scale[set] += rate_scale_cat_no_y_branch_id(c,b)*p;
      }
    }

    static const smlfloat expect_tolerance(1.e-6);
    for(unsigned set(0);set<n_sets;++set){
      expect_scale[set] /= set_prob[set];
      if(std::fabs(expect_scale[set]-expect)>expect_tolerance){
        std::ostringstream oss;
        oss << "mutation model set normalization failed: "
            << "expect,n_sets,set,expect_scale[set],prob[set]: "
            << expect << " " << n_sets << " " << set << " "
            << expect_scale[set] << " " << set_prob[set] << "\n";
        throw substk_exception(oss.str().c_str());
      }
    }
  }

#if 0
  if(is_time_scaled){
    const unsigned ntc(typed_cat_size(CAT_PARAM_TYPE::TIME));
    for(unsigned i(0);i<ntc;++i){
     const unsigned set(cm.mut_model_norm_set_from_time_cat(i));
     time_cat_time_scale(i,1./mut_model_scale[set]);
    }
  }
#endif
}



void
rate_gtor_nuc_base::
fix_input_param_internal(){

  base_t::fix_input_param_internal();

  {
    const CAT_PARAM_TYPE::index_t pt(CAT_PARAM_TYPE::MUT_RATE);
    CAT_MIX_TYPE::index_t mt(CAT_MIX_TYPE::SITE);
    fix_input_param_cat_branch_paramset(PARAMSET_SITE_RATE_CATS,pt,mt,
                                        cm().mut_rate_norm_set_size(mt),
                                        cm().mut_rate_norm_set_from_mut_rate_cat(mt));
    mt=CAT_MIX_TYPE::GROUP;
    fix_input_param_cat_branch_paramset(PARAMSET_GROUP_RATE_CATS,pt,mt,
                                        cm().mut_rate_norm_set_size(mt),
                                        cm().mut_rate_norm_set_from_mut_rate_cat(mt));
  }
  normalize_mut_rate_params();
}



/// \todo this version needs to scale the tree times somehow??
void
rate_gtor_nuc_base::
param_update_from_new_bg_pdistro_internal(){
  base_t::param_update_from_new_bg_pdistro_internal();
  normalize_mut_rate_params();
}



unsigned
rate_gtor_nuc_base::
paramset_size(const int p) const {
  switch(static_cast<const paramset_t>(p)){
  case PARAMSET_SITE_RATE_CATS :
    return cm().typed_cat_size(CAT_PARAM_TYPE::MUT_RATE,CAT_MIX_TYPE::SITE);
  case PARAMSET_GROUP_RATE_CATS :
    return cm().typed_cat_size(CAT_PARAM_TYPE::MUT_RATE,CAT_MIX_TYPE::GROUP);
  case PARAMSET_NUC :
    {
      unsigned nsize(0);
      const unsigned nmmc(cm().typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL));
      for(unsigned i(0);i<nmmc;++i){
        nsize += paramset_mutation_model_cat_size(i);
      }
      return nsize;
    }
  case PARAMSET_EDGE_STRENGTH :
    if(is_use_edge_correction()) {return 2;}
    else                         {return 0;}
  default : return 0;
  }
}




void
rate_gtor_nuc_base::
paramset_label(const paramset_t p,
               const unsigned n,
               std::ostream& os) const {

  switch(p){
  case PARAMSET_SITE_RATE_CATS :
    os << "site_rate_cats_" << n;
    break;
#if 0
  case PARAMSET_SITE_RATE_CATS_PROB :
    os << "site_rate_cats_prior_" << n;
    break;
  case PARAMSET_SITE_MUTATION_MODEL_CATS_PROB :
    os << "site_mutation_model_prior_" << n;
    break;
#endif
  case PARAMSET_GROUP_RATE_CATS :
    os << "group_rate_cats_" << n;
    break;
#if 0
  case PARAMSET_GROUP_RATE_CATS_PROB :
    os << "group_rate_cats_prior_" << n;
    break;
#endif
  case PARAMSET_NUC :
    {
      bool is_valid(false);
      unsigned mut_model_no(0);
      unsigned model_param_no(0);
      {
        const unsigned nmmc(cm().typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL));
        int ncopy(n);
        for(unsigned i(0);i<nmmc;++i){
          const int s(paramset_mutation_model_cat_size(i));
          if(ncopy<s){
            if(ncopy>=0){
              mut_model_no=i;
              model_param_no=ncopy;
              is_valid=true;
            }
            break;
          }
          ncopy -= s;
        }
      }
      if(! is_valid) die("invalid call to paramset_label");
      get_nuc_param_index_name(mut_model_no,model_param_no,os);
    }
    break;
  case PARAMSET_EDGE_STRENGTH :
    os << "edge_strength_" << n;
    break;
  default :
    os << "unknown_"  << n;
  }
}



const char* const section_id = "rate_gtor_nuc_base";
const char* const end_label = "END";

const iliner il(section_id);


void
rate_gtor_nuc_base::
store_state(std::ostream& os) const {

  base_t::store_state(os);

  os << section_id << " is_use_edge_correction " << is_use_edge_correction() << "\n";
  const unsigned nmmc(cm().typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL));
  for(unsigned i(0);i<nmmc;++i){
    os << section_id << " rate_model_" << i << " " << rate_model_nuc(i) << "\n";
  }
  for(unsigned i(0);i<nmmc;++i){
    os << section_id << " context_model_" << i << " " << context_model_nuc(i) << "\n";
  }

  // write all paramsets:
  //
  for(unsigned p(0);p<PARAMSET_END;++p){
    const unsigned st=paramset_start(static_cast<paramset_t>(p));
    const unsigned s=paramset_size(static_cast<paramset_t>(p));
    for(unsigned i(0);i<s;++i){
      os << section_id << " ";
      paramset_label(static_cast<paramset_t>(p),i,os);
      os << " " << _is_train_param[st+i] << " " << _param[st+i] << "\n";
    }
  }

  os << section_id << " " << end_label << "\n";
}



void
rate_gtor_nuc_base::
load_state_internal(std::istream& is) {

  base_t::load_state_internal(is);

  il.advance(is,"is_use_edge_correction");
  is >> _data.nopt.is_use_edge_correction;

  const unsigned nmmc(cm().typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL));
  for(unsigned i(0);i<nmmc;++i){
    std::ostringstream oss;
    oss << "rate_model_" << i;
    il.advance(is,oss.str().c_str());
    const RATE_MODEL_NUC::index_t rm(static_cast<RATE_MODEL_NUC::index_t>(read_int(is)));
    _data.nopt.set_rate_model(i,rm);
  }

  for(unsigned i(0);i<nmmc;++i){
    std::ostringstream oss;
    oss << "context_model_" << i;
    il.advance(is,oss.str().c_str());
    const CONTEXT_MODEL_NUC::index_t cmn(static_cast<CONTEXT_MODEL_NUC::index_t>(read_int(is)));
    _data.nopt.set_context_model(i,cmn);
  }

  init_class_storage();

  // read all paramsets:
  //
  for(unsigned p(0);p<PARAMSET_END;++p){
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




// helper func to set_ts_tv_ratio
void
rate_gtor_nuc_base::
set_ts_tv_ratio_nhood(const unsigned mut_model_cat,
                      const smlfloat ratio,
                      const NUC::index_t n0,
                      const NUC::index_t n2){
  for(unsigned n1(0);n1<NUC::SIZE;++n1){
    for(unsigned to(0);to<NUC::SIZE;++to){
      if(n1==to) continue;

      const NUC::index_t n1x = static_cast<NUC::index_t>(n1);
      const NUC::index_t tox = static_cast<NUC::index_t>(to);
      if(context_model_nuc(mut_model_cat) == CONTEXT_MODEL_NUC::FACTORED_TRIPLET){
        const int i1(nuc_context_mut_rate_param_index(n0,n1x,n2,tox,mut_model_cat,JOINT1));
        if(i1>=0) _param[i1] = ( NUC::is_transition(n1x,tox) ? std::sqrt(ratio) : 1. );

        const int i2(nuc_context_mut_rate_param_index(n0,n1x,n2,tox,mut_model_cat,JOINT2));
        if(i2>=0) _param[i2] = ( NUC::is_transition(n1x,tox) ? std::sqrt(ratio) : 1. );

      } else {
        const int i(nuc_context_mut_rate_param_index(n0,n1x,n2,tox,mut_model_cat));
        if(i>=0) _param[i] = ( NUC::is_transition(n1x,tox) ? ratio : 1. );
      }
    }
  }
}


void
rate_gtor_nuc_base::
set_ts_tv_ratio(const unsigned mut_model_cat,
                const smlfloat ratio){

  const unsigned pst(paramset_mutation_model_cat_start(mut_model_cat));
  const unsigned ps(paramset_mutation_model_cat_size(mut_model_cat));
  for(unsigned i(pst);i<(pst+ps);++i)  _param[i] = 1;

  using namespace CONTEXT_GROUP;
  const index_t cg(context_model_group(context_model_nuc(mut_model_cat)));

  if       (cg == NONE){
    set_ts_tv_ratio_nhood(mut_model_cat,ratio,NUC::N,NUC::N);
  } else if(cg == ONLY_PRE){
    for(unsigned n0(0);n0<NUC::SIZE;++n0){
      set_ts_tv_ratio_nhood(mut_model_cat,ratio,
                            static_cast<NUC::index_t>(n0),
                            NUC::N);
    }
  } else if(cg == ONLY_POST){
    for(unsigned n2(0);n2<NUC::SIZE;++n2){
      set_ts_tv_ratio_nhood(mut_model_cat,ratio,
                            NUC::N,
                            static_cast<NUC::index_t>(n2));
    }
  } else if(cg == TWOWAY){
    for(unsigned n0(0);n0<NUC::SIZE;++n0){
      for(unsigned n2(0);n2<NUC::SIZE;++n2){
        set_ts_tv_ratio_nhood(mut_model_cat,ratio,
                              static_cast<NUC::index_t>(n0),
                              static_cast<NUC::index_t>(n2));
      }
    }
  } else { die("Unsupported nuc context model"); }
}


/// non-reversible param coding:
///    A  C  G  T
/// A  -  0  1  2
/// C  3  -  4  5
/// G  6  7  -  8
/// T  9 10 11  -
///
/// reversible param coding:
///    A  C  G  T
/// A  -  0  1  2
/// C  0  -  3  4
/// G  1  3  -  5
/// T  2  4  5  -
///


//          i2
//        0 1 2
//      +------
//    0 | - 0 1
// i1 1 | - - 2
//    2 | - - -
//
static
unsigned
symm_matrix_nodiag_index(unsigned i1,
                             unsigned i2,
                             const unsigned size){

  assert(i1!=i2 && i1<size && i2<size);
  if(i2<i1) std::swap(i2,i1);
  return (i2-1-i1)+(size-1)*i1-(i1*(i1-1))/2;
}
#if 0
static
void
symm_matrix_nodiag_index_inverse(const unsigned i,
                                 const unsigned size,
                                 unsigned& i1,
                                 unsigned& i2){

  assert(i<(off_diag_size(size)/2));
  abort(); /// \todo finish this
}
#endif



//          i2
//        0 1 2
//      +------
//    0 | - 0 1
// i1 1 | 2 - 3
//    2 | 4 5 -
//
static
unsigned
square_matrix_nodiag_index(const unsigned i1,
                           const unsigned i2,
                           const unsigned size){

  assert(i1!=i2 && i1<size && i2<size);
  return i2+(i2<i1?0:-1)+(size-1)*i1;
}

static
void
square_matrix_nodiag_index_inverse(const unsigned i,
                                   const unsigned size,
                                   unsigned& i1,
                                   unsigned& i2){

  assert(i<off_diag_size(size));
  i1=i/(size-1);
  i2=i%(size-1);
  if(i1<=i2) i2 += 1;
}



static
void
get_nuc_param_index_name_nocontext(const RATE_MODEL_NUC::index_t rmn,
                                   const unsigned param_index,
                                   std::ostream& os) {

  {
    using namespace RATE_MODEL_NUC;

    assert(param_index<get_rate_model_param_size(rmn));

    if       (rmn == JC69 || rmn == F81){
      os << "nuc_mut";
    } else if(rmn == K80 || rmn == HKY85 || rmn == GY94){
      if(param_index==1){ os << "transition"; }
      else              { os << "transversion"; }
    } else if(rmn == REV){
      os << "X1<->X2";
    } else if(rmn == NONREV){
      unsigned n1,n2;
      square_matrix_nodiag_index_inverse(param_index,NUC::SIZE,n1,n2);
      os << NUC::syms[n1] << "->" << NUC::syms[n2];
    } else {
      die("Unsupported nuc rate model");
    }
    os << "_rate";
  }
}




// produce a human readable name for the nucleotide parameters
//
void
rate_gtor_nuc_base::
get_nuc_param_index_name(const unsigned mut_model_cat_no,
                         const unsigned param_index,
                         std::ostream& os) const {

  os << "nuc_mutation_model_" << mut_model_cat_no << "_param_";

  const CONTEXT_MODEL_NUC::index_t cmn(context_model_nuc(mut_model_cat_no));
  const RATE_MODEL_NUC::index_t rmn(rate_model_nuc(mut_model_cat_no));

  assert(param_index<get_rate_context_model_param_size(rmn,cmn));

  {
    static const char cpg_tag[] = "CpG";

    using namespace CONTEXT_MODEL_NUC;

    if       (cmn == INDY){
      get_nuc_param_index_name_nocontext(rmn,param_index,os);
      os << "_indy_context";
    } else if(cmn == PRE_DOUBLET){
      get_nuc_param_index_name_nocontext(rmn,param_index/NUC::SIZE,os);
      const unsigned from0(param_index%NUC::SIZE);
      os << "_" << NUC::syms[from0] << "N_context";
    } else if(cmn == POST_DOUBLET){
      get_nuc_param_index_name_nocontext(rmn,param_index/NUC::SIZE,os);
      const unsigned from2(param_index%NUC::SIZE);
      os << "_N" << NUC::syms[from2] << "_context";
    } else if(cmn == FACTORED_TRIPLET){
      get_nuc_param_index_name_nocontext(rmn,param_index/(NUC::SIZE+NUC::SIZE),os);
      os << "_xx_context";
    } else if(cmn == TRIPLET){
      get_nuc_param_index_name_nocontext(rmn,param_index/(NUC::SIZE*NUC::SIZE),os);
      const unsigned tmpu(param_index%(NUC::SIZE*NUC::SIZE));
      const unsigned from0(tmpu%NUC::SIZE);
      const unsigned from2(tmpu/NUC::SIZE);
      os << "_" << NUC::syms[from0] << "N" << NUC::syms[from2] <<"_context";
    } else if(cmn == CPG_ONLY || cmn == TPA_ONLY){
      static const char tpa_tag[] = "TpA";

      const char* mut_tag(cmn == CPG_ONLY ? cpg_tag : tpa_tag);

      const unsigned psn(paramset_mutation_model_cat_size(mut_model_cat_no));
      if     (param_index == psn-2){ os << mut_tag << "_transversion_rate"; }
      else if(param_index == psn-1){ os << mut_tag << "_transition_rate"; }
      else { get_nuc_param_index_name_nocontext(rmn,param_index,os); }
    } else if(cmn == CPG_1TI){

      const unsigned psn(paramset_mutation_model_cat_size(mut_model_cat_no));
      if(param_index == psn-1){ os << cpg_tag << "_transition_rate"; }
      else { get_nuc_param_index_name_nocontext(rmn,param_index,os); }
    } else if(cmn == CPG_2TI){

      const unsigned psn(paramset_mutation_model_cat_size(mut_model_cat_no));
      if     (param_index == psn-2){ os << cpg_tag << "_C->T_transition_rate"; }
      else if(param_index == psn-1){ os << cpg_tag << "_G->A_transition_rate"; }
      else { get_nuc_param_index_name_nocontext(rmn,param_index,os); }
    } else if(cmn == CPG_NONREV){
      const unsigned psn(paramset_mutation_model_cat_size(mut_model_cat_no));
      os << cpg_tag;
      if     (param_index == psn-1){ os << "_G->A_rate"; }
      else if(param_index == psn-2){ os << "_G->C_rate"; }
      else if(param_index == psn-3){ os << "_G->T_rate"; }
      else if(param_index == psn-4){ os << "_C->A_rate"; }
      else if(param_index == psn-5){ os << "_C->G_rate"; }
      else if(param_index == psn-6){ os << "_C->T_rate"; }

      else { get_nuc_param_index_name_nocontext(rmn,param_index,os); }



    } else { die("Unsupported nuc context model"); }
  }

  os << "_" << param_index;
}




static
unsigned
get_rate_model_nuc_mut_rate_param_index(const RATE_MODEL_NUC::index_t i,
                                        const NUC::index_t from,
                                        const NUC::index_t to){

  using namespace RATE_MODEL_NUC;

  switch(i){
  case JC69:
  case F81:    return 0;
  case K80:
  case HKY85:
  case GY94:
    if( NUC::is_transition(from,to) ) return 1;
    else                              return 0;
  case REV:    return symm_matrix_nodiag_index(from,to,NUC::SIZE);
  case NONREV: return square_matrix_nodiag_index(from,to,NUC::SIZE);
  default: die("Unsupported nuc rate model");
  }
}




// designed for context-dependent case, positions
// are numbered 012 and nucleotide change occurs
// at position 1
//
unsigned
rate_gtor_nuc_base::
nuc_context_mut_rate_param_index(const NUC::index_t from0,
                                 const NUC::index_t from1,
                                 const NUC::index_t from2,
                                 const NUC::index_t to,
                                 const unsigned mut_model_cat,
                                 const joint_param_t jp) const {

  if( from1 == to || from1 == NUC::N || to == NUC::N) {
    die("Invalid parameters for rate model");
  }

  const unsigned psn(paramset_mutation_model_cat_size(mut_model_cat));
  unsigned param_index(get_rate_model_nuc_mut_rate_param_index(rate_model_nuc(mut_model_cat),from1,to));

  { // modify param_index according to context:

    using namespace CONTEXT_MODEL_NUC;
    const index_t i(context_model_nuc(mut_model_cat));

    if((from0 == NUC::N) &&
       ((i == PRE_DOUBLET) ||
        (i == FACTORED_TRIPLET && jp == JOINT1) ||
        (i == TRIPLET))) die("get_nuc_mut_rate_param_index: invalid from0 nuc");

    if((from2 == NUC::N) &&
       ((i == POST_DOUBLET) ||
        (i == FACTORED_TRIPLET && jp == JOINT2) ||
        (i == TRIPLET))) die("get_nuc_mut_rate_param_index: invalid from2 nuc");

    if       (i == INDY){
      // pass
    } else if(i == PRE_DOUBLET){
      param_index = param_index*(NUC::SIZE)+from0;
    } else if(i == POST_DOUBLET){
      param_index = param_index*(NUC::SIZE)+from2;
    } else if(i == FACTORED_TRIPLET){
      if(jp == JOINT1) { param_index = param_index*(NUC::SIZE+NUC::SIZE)+from0; }
      else             { param_index = param_index*(NUC::SIZE+NUC::SIZE)+(NUC::SIZE+from2); }
    } else if(i == TRIPLET){
      param_index = param_index*(NUC::SIZE*NUC::SIZE)+(from2*NUC::SIZE+from0);
    } else if(i == CPG_ONLY){
      if((from0 == NUC::C && from1 == NUC::G) ||
         (from1 == NUC::C && from2 == NUC::G)) {
        if (NUC::is_transition(from1,to)) { param_index = psn-1; }
        else                              { param_index = psn-2; }
      }
    } else if(i == TPA_ONLY){
      if((from0 == NUC::T && from1 == NUC::A) ||
         (from1 == NUC::T && from2 == NUC::A)) {
        if (NUC::is_transition(from1,to)) { param_index = psn-1; }
        else                              { param_index = psn-2; }
      }
    } else if(i == CPG_1TI) {
      if((from0 == NUC::C && from1 == NUC::G) ||
         (from1 == NUC::C && from2 == NUC::G)) {
        if(NUC::is_transition(from1,to)){ param_index = psn-1; }
      }
    } else if(i == CPG_2TI) {
      if(NUC::is_transition(from1,to)){
        if     (from0 == NUC::C && from1 == NUC::G){ param_index = psn-1; }
        else if(from1 == NUC::C && from2 == NUC::G){ param_index = psn-2; }
      }
    } else if(i == CPG_NONREV) {
        if       (from0 == NUC::C && from1 == NUC::G){
          if     (to == NUC::A){ param_index = psn-1; }
          else if(to == NUC::C){ param_index = psn-2; }
          else if(to == NUC::T){ param_index = psn-3; }
        } else if(from1 == NUC::C && from2 == NUC::G){
          if     (to == NUC::A){ param_index = psn-4; }
          else if(to == NUC::G){ param_index = psn-5; }
          else if(to == NUC::T){ param_index = psn-6; }
        }

    } else { die("Unsupported nuc context model"); }
  }

  assert(param_index < psn);
  return param_index+paramset_mutation_model_cat_start(mut_model_cat);
}



// get nucleotide parameter from NUC codes... this function layer is
// mainly here to handle the factored parameter case:
//
irv_t<smlfloat>
rate_gtor_nuc_base::
nuc_context_mut_rate_param(const NUC::index_t from0,
                           const NUC::index_t from1,
                           const NUC::index_t from2,
                           const NUC::index_t to,
                           const unsigned mut_model_cat,
                           const trainable_mode::index_t tm) const {

  const CONTEXT_GROUP::index_t cg(context_model_group(context_model_nuc(mut_model_cat)));

  if(((cg == CONTEXT_GROUP::ONLY_PRE  || cg == CONTEXT_GROUP::TWOWAY) && (from0 == NUC::N)) ||
     ((cg == CONTEXT_GROUP::ONLY_POST || cg == CONTEXT_GROUP::TWOWAY) && (from2 == NUC::N))){
    die("function requires known context");
  }

  if(context_model_nuc(mut_model_cat) == CONTEXT_MODEL_NUC::FACTORED_TRIPLET){
    const unsigned param_index1(nuc_context_mut_rate_param_index(from0,from1,from2,to,mut_model_cat,JOINT1));
    const irv_t<smlfloat> param1(_param[param_index1],param_variance(param_index1));

    const unsigned param_index2(nuc_context_mut_rate_param_index(from0,from1,from2,to,mut_model_cat,JOINT2));
    const irv_t<smlfloat> param2(_param[param_index2],param_variance(param_index2));

    const bool is_full_joint(_is_train_param[param_index1] && _is_train_param[param_index2]);
    const bool is_half_joint((!is_full_joint) &&
                             (_is_train_param[param_index1] || _is_train_param[param_index2]));

    const bool is_skip((tm == trainable_mode::HALF_JOINT && (! is_half_joint)) ||
                       (tm == trainable_mode::FULL && (! is_full_joint)));

    if(is_skip){
      return irv_t<smlfloat>(0.);
    } else {
      return param1*param2;
    }

  } else {
    const unsigned param_index(nuc_context_mut_rate_param_index(from0,from1,from2,to,mut_model_cat));
    const bool is_skip((tm == trainable_mode::HALF_JOINT) ||
                       (tm == trainable_mode::FULL && (! _is_train_param[param_index])));

    if(is_skip){
      return irv_t<smlfloat>(0.);
    } else {
      return irv_t<smlfloat>(_param[param_index],param_variance(param_index));
    }
  }
}



///
//
irv_t<smlfloat>
rate_gtor_nuc_base::
nuc_context_mut_rate(const NUC::index_t from0,
                     const NUC::index_t from1,
                     const NUC::index_t from2,
                     const NUC::index_t to,
                     const unsigned cat_no,
                     const unsigned branch_cat_set,
                     const prob_t* nuc_bg0,
                     const prob_t* nuc_bg2,
                     const trainable_mode::index_t tm) const {

  const unsigned mut_model_cat(cm().typed_cat_no_from_cat_no_y_branch_cat_set(cat_no,branch_cat_set,CAT_PARAM_TYPE::MUT_MODEL));
  const CONTEXT_GROUP::index_t cg(context_model_group(context_model_nuc(mut_model_cat)));

  if       ((cg == CONTEXT_GROUP::ONLY_PRE || cg == CONTEXT_GROUP::TWOWAY) && (from0 == NUC::N)){
    if(nuc_bg0==0) die("unspecified edge nuc distro");
    irv_t<smlfloat> param(0);
    for(unsigned i(0);i<NUC::SIZE;++i){
      const irv_t<smlfloat> ptmp = nuc_context_mut_rate(static_cast<NUC::index_t>(i),from1,from2,to,
                                                        cat_no,branch_cat_set,nuc_bg0,nuc_bg2,tm);
      param += ptmp*nuc_bg0[i];
    }
    return param;

  } else if((cg == CONTEXT_GROUP::ONLY_POST || cg == CONTEXT_GROUP::TWOWAY) && (from2 == NUC::N)){
    if(nuc_bg2==0) die("unspecified edge nuc distro");
    irv_t<smlfloat> param(0);
    for(unsigned i(0);i<NUC::SIZE;++i){
      const irv_t<smlfloat> ptmp = nuc_context_mut_rate(from0,from1,static_cast<NUC::index_t>(i),to,
                                                        cat_no,branch_cat_set,nuc_bg0,nuc_bg2,tm);
      param += ptmp*nuc_bg2[i];
    }
    return param;

  } else {
    irv_t<smlfloat> param(nuc_context_mut_rate_param(from0,from1,from2,to,mut_model_cat,tm));

    {
      using namespace RATE_MODEL_NUC;
      const index_t rmn(rate_model_nuc(mut_model_cat));
      if(rmn == F81 || rmn == HKY85 || rmn == REV){
        return param*bg_pdistro_nuc_cat(cat_no)[to];
      } else {
        return param;
      }
    }
  }
}



smlfloat
rate_gtor_nuc_base::
cat_nuc_context_mut_rate(const NUC::index_t from0,
                         const NUC::index_t from1,
                         const NUC::index_t from2,
                         const NUC::index_t to,
                         const unsigned cat_no,
                         const unsigned branch_cat_set,
                         const prob_t* nuc_bg0,
                         const prob_t* nuc_bg2,
                         const bool is_coding_model) const {

  const unsigned mmc(cm().typed_cat_no_from_cat_no_y_branch_cat_set(cat_no,branch_cat_set,CAT_PARAM_TYPE::MUT_MODEL));
  if(rate_model_nuc(mmc)==RATE_MODEL_NUC::GY94 && ! is_coding_model){
    die("GY94 rates are valid for coding models only");
  }

  smlfloat val(nuc_context_mut_rate(from0,from1,from2,to,cat_no,branch_cat_set,nuc_bg0,nuc_bg2));

  val *= get_paramset_member(PARAMSET_SITE_RATE_CATS,
                             cm().typed_cat_no_from_cat_no_y_branch_cat_set(cat_no,branch_cat_set,CAT_PARAM_TYPE::MUT_RATE,CAT_MIX_TYPE::SITE));
  val *= get_paramset_member(PARAMSET_GROUP_RATE_CATS,
                             cm().typed_cat_no_from_cat_no_y_branch_cat_set(cat_no,branch_cat_set,CAT_PARAM_TYPE::MUT_RATE,CAT_MIX_TYPE::GROUP));

  return val;
}



/// \todo eventually, this method should never touch locked parameters,
/// they should be initialized together with init_param_locks
///
void
rate_gtor_nuc_base::
reset_param_internal(const PARAM_INIT_TYPE::index_t pinit) {

  base_t::reset_param_internal(pinit);

  if(pinit == PARAM_INIT_TYPE::RANDOM){
    static const smlfloat rate_cats_min(0.1);
    static const smlfloat rate_cats_range(2.);

    {
      const unsigned st(paramset_start(PARAMSET_SITE_RATE_CATS));
      const unsigned s(paramset_size(PARAMSET_SITE_RATE_CATS));
      for(unsigned i(st);i<(st+s);++i){
        if(_is_train_param[i]) { _param[i] = random_uniform()*rate_cats_range+rate_cats_min; }
        else                   { _param[i] = 1.; }
      }
    }
    {
      const unsigned st(paramset_start(PARAMSET_GROUP_RATE_CATS));
      const unsigned s(paramset_size(PARAMSET_GROUP_RATE_CATS));
      for(unsigned i(st);i<(st+s);++i){
        if(_is_train_param[i]) { _param[i] = random_uniform()*rate_cats_range+rate_cats_min; }
        else                   { _param[i] = 1.; }
      }
    }

  } else if(pinit == PARAM_INIT_TYPE::START){
    static const smlfloat start_ratio(3.);

    const unsigned nmmc(cm().typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL));
    for(unsigned i(0);i<nmmc;++i){
      set_ts_tv_ratio(i,start_ratio);
    }

    // don't setup initial category parameters on a slope if they are
    // tied to assigned data
    if(cm().assigned_data_set_size() > 1){ die("not yet setup to init parameters for multi adsets"); }
    set_paramset_slope(PARAMSET_GROUP_RATE_CATS);
    set_paramset_slope(PARAMSET_SITE_RATE_CATS);
    set_mutation_model_slope();
  }

  {  // for any init case...
    const unsigned st(paramset_start(PARAMSET_EDGE_STRENGTH));
    const unsigned s(paramset_size(PARAMSET_EDGE_STRENGTH));
    for(unsigned i(st);i<(st+s);++i) if(_is_train_param[i]) _param[i]=0.5;
  }
}



static
void
bg_pdistro_nuc_conditioned_update(const SITE_MODEL::index_t sm,
                                  const prob_t* bg_pdistro,
                                  prob_t* bg_pdistro_nuc_conditioned_on_3p,
                                  prob_t* bg_pdistro_nuc_conditioned_on_5p){

  get_nuc_distro_conditioned_on_3p_from_site_model_distro(sm,
                                                          bg_pdistro,
                                                          bg_pdistro_nuc_conditioned_on_3p);

  get_nuc_distro_conditioned_on_5p_from_site_model_distro(sm,
                                                          bg_pdistro,
                                                          bg_pdistro_nuc_conditioned_on_5p);
}



void
rate_gtor_nuc_base::
bg_pdistro_nuc_conditioned_update(){
  const SITE_MODEL::index_t sm(site_model());

  const unsigned bcn(bg_cat_size());
  _data.bg_pdistro_nuc_conditioned_on_3p_bg_cat.init(bcn,NUC::SIZE*NUC::SIZE);
  _data.bg_pdistro_nuc_conditioned_on_5p_bg_cat.init(bcn,NUC::SIZE*NUC::SIZE);
  for(unsigned c(0);c<bcn;++c){
    ::bg_pdistro_nuc_conditioned_update(sm,bg_pdistro_bg_cat(c),
                                        _data.bg_pdistro_nuc_conditioned_on_3p_bg_cat[c],
                                        _data.bg_pdistro_nuc_conditioned_on_5p_bg_cat[c]);
  }
}



/// helper func to report():
void
rate_gtor_nuc_base::
report_rate_set(const NUC::index_t n0,
                const NUC::index_t n2,
                const unsigned cat_no,
                const unsigned branch_cat_set,
                std::ostream& os) const {

  const char sn0(NUC::syms[n0]);
  const char sn2(NUC::syms[n2]);

  const unsigned mmc(cm().typed_cat_no_from_cat_no_y_branch_cat_set(cat_no,branch_cat_set,CAT_PARAM_TYPE::MUT_MODEL));

  std::string s("MUT RATE");
  if(context_model_nuc(mmc) != CONTEXT_MODEL_NUC::INDY){
    s = s+" : context= "+sn0+"X"+sn2+"->Z\n";
  }
  irv_t<smlfloat> rts[NUC::SIZE*NUC::SIZE];
  for(unsigned i(0);i<NUC::SIZE;++i){
    const NUC::index_t n1(static_cast<NUC::index_t>(i));
    for(unsigned j(0);j<NUC::SIZE;++j){
      if(i==j) {
        rts[j+i*NUC::SIZE] = irv_t<smlfloat>(0.);
        continue;
      }
      const NUC::index_t n1x(static_cast<NUC::index_t>(j));

      rts[j+i*NUC::SIZE] = cat_nuc_context_mut_rate(n0, n1,
                                                    n2, n1x,
                                                    cat_no,
                                                    branch_cat_set,
                                                    bg_pdistro_nuc_conditioned_on_3p_cat(cat_no,n1),
                                                    bg_pdistro_nuc_conditioned_on_5p_cat(cat_no,n1));
    }
  }

  os << s << "\n";

  unsigned cell_width(9);
  if(rts[1].v != 0.){ cell_width = 22; }
  static const unsigned prec(3);

  matrix_report(rts,NUC::SIZE,NUC::syms,os,cell_width,prec);


  // print out a higher precision list in easily parsed form:
  static const unsigned lprec(6);
  os << std::setprecision(lprec);
  for(unsigned i(0);i<NUC::SIZE;++i){
    const char sni(NUC::syms[i]);
    for(unsigned j(0);j<NUC::SIZE;++j){
      if(i==j) continue;
      const char snj(NUC::syms[j]);
      os << "CDM_RATE: " <<  sn0 << sni << sn2 << "->" << snj << " " << rts[j+i*NUC::SIZE] << "\n";
    }
  }

  // get transi/tranv:
  //
  {
    irv_t<smlfloat> transi,transv;
    get_si_sv(transi,transv,rts,bg_pdistro_nuc_cat(cat_no));
    os << "CT_TI: " << transi << "\n";
    os << "CT_TV: " << transv << "\n";
    os << "CT_TI/TV: " << transi/transv << "\n\n";
  }

  // get W->S/S->W
  //
  const irv_t<smlfloat> wsr(get_ws_ratio(rts,bg_pdistro_nuc_cat(cat_no)));
  os << "CT_WS/SW: " << wsr << "\n";
}



static
void
print_section_divider(std::ostream& os){
  os << "######################################################################\n";
}



static
void
print_branch_cat_mut_model_header(const unsigned branch_set_id,
                                  const std::string& branch_set_label,
                                  const unsigned cat_id,
                                  const std::string& cat_label,
                                  const unsigned mut_model_id,
                                  const std::string& mut_model_label,
                                  std::ostream& os){

  print_section_divider(os);
  os << "branch_cat_set_seq_cat_mut_model_summary:\n\n";
  os << "branch_cat_set: " << branch_set_id << " " << branch_set_label << "\n";
  os << "global_seq_cat: " << cat_id << " " << cat_label << "\n";
  os << "mut_model_param_cat: " << mut_model_id << " " << mut_model_label << "\n\n";
}



void
rate_gtor_nuc_base::
report_rate_cats(const CAT_MIX_TYPE::index_t mt,
                 const int paramset,
                 const unsigned branch_id,
                 std::ostream& os) const {

  const unsigned nrc(cm().typed_cat_size(CAT_PARAM_TYPE::MUT_RATE,mt));

  simple_array<prob_t> rc_prior(nrc);
  cm().branch_typed_cat_pdistro(rc_prior.ptr(),branch_id,CAT_PARAM_TYPE::MUT_RATE,mt);

  std::vector<std::pair<irv_t<smlfloat>,std::pair<irv_t<smlfloat>,std::string> > > tmp_rate_cats(nrc);

  for(unsigned i(0);i<nrc;++i) {
    const std::string& cat_label(cm().typed_cat_label(i,CAT_PARAM_TYPE::MUT_RATE,mt));
    tmp_rate_cats[i] =
      std::make_pair(get_paramset_member(paramset,i),
                     std::make_pair(rc_prior[i],cat_label));
  }

  std::sort(tmp_rate_cats.begin(),tmp_rate_cats.end());

  os << std::string(CAT_MIX_TYPE::syms[mt]) << " rate categories (" << nrc << "):\n";
  os.unsetf(std::ios::fixed);
  os << std::setprecision(7);
  for(unsigned i(0);i<nrc;++i){
    os << "[ " << tmp_rate_cats[i].second.second << " ]  :"
       << " scale=" << tmp_rate_cats[i].first
       << " branch_prior_prob=(" << tmp_rate_cats[i].second.first << ")\n";
  }
  os << "\n";
}



void
rate_gtor_nuc_base::
report(std::ostream& os) const {

  using namespace CONTEXT_MODEL_NUC;

  const unsigned n_cats(cm().cat_size());

  for(unsigned c(0);c<n_cats;++c){
    const std::string& cat_label(cm().cat_label(c));
    const rates_func_options_base bopt(c);

    const unsigned n_branch_cat_sets(cm().branch_cat_set_size(c));
    for(unsigned bcs(0);bcs<n_branch_cat_sets;++bcs){
      std::string branch_set_label;
      cm().branch_cat_set_label(c,bcs,branch_set_label);

      os << "cat: " << cat_label << " branch_cat_set: " << branch_set_label << "\n";

      const rates_func_options opt(bopt,bcs);

      smlfloat stat_pdistro[NUC::SIZE];
      stat_pdistro_nuc(stat_pdistro,opt);
      os << "stat_nuc_distro:\n";
      for(unsigned i(0);i<NUC::SIZE;++i){
        os << NUC::syms[i] << ": " << stat_pdistro[i] << "\n";
      }
      os << "bg_nuc_distro:\n";
      for(unsigned i(0);i<NUC::SIZE;++i){
        os << NUC::syms[i] << ": " << bg_pdistro_nuc_cat(c)[i] << "\n";
      }
      os << "\n\n";
    }
  }

  simple_array<prob_t> cp(n_cats);
  cm().cat_pdistro(cp.ptr());

  const unsigned n_branches(cm().branch_size());
  for(unsigned b(0);b<n_branches;++b){
    print_section_divider(os);
    os << "Reporting mutation cat summary for branch [" << cm().branch_label(b) << "]\n\n";

    {  // mutation model cats
      const unsigned nmmc(cm().typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL));

      simple_array<prob_t> mmc_prior(nmmc);
      cm().branch_typed_cat_pdistro(mmc_prior.ptr(),b,CAT_PARAM_TYPE::MUT_MODEL);

      std::vector<std::pair<irv_t<smlfloat>,std::pair<irv_t<smlfloat>,std::string> > > tmp_rate_cats(nmmc);
      for(unsigned i(0);i<nmmc;++i) {
        smlfloat scale_val(0.);
        prob_t pi(0.);
        for(unsigned c(0);c<n_cats;++c){
          const unsigned mmc(cm().typed_cat_no_from_cat_no_y_branch_id(c,b,CAT_PARAM_TYPE::MUT_MODEL));
          if(mmc != i) continue;
          scale_val += rate_scale_cat_no_y_branch_id(c,b)*cp[c];
          pi += cp[c];
        }

        scale_val /= pi;

        const std::string& cat_label(cm().typed_cat_label(i,CAT_PARAM_TYPE::MUT_MODEL));
        tmp_rate_cats[i] =
          std::make_pair(scale_val,
                         std::make_pair(mmc_prior[i],cat_label));
      }

      std::sort(tmp_rate_cats.begin(),tmp_rate_cats.end());

      os << "mutation model categories (" << nmmc << "):\n";
      os.unsetf(std::ios::fixed);
      os << std::setprecision(7);
      for(unsigned i(0);i<nmmc;++i){
        os << "[ " << tmp_rate_cats[i].second.second << " ]  :"
           << " scale=" << tmp_rate_cats[i].first
           << " branch_prior_prob=(" << tmp_rate_cats[i].second.first << ")\n";
      }
      os << "\n";
    }

    report_rate_cats(CAT_MIX_TYPE::SITE,PARAMSET_SITE_RATE_CATS,b,os);

    report_rate_cats(CAT_MIX_TYPE::GROUP,PARAMSET_GROUP_RATE_CATS,b,os);
  }


  /// \todo make this report less verbose for cases when cats are distinguished by
  ///       selection parameters only
  ///
  for(unsigned c(0);c<n_cats;++c){
    const std::string& cat_label(cm().cat_label(c));

    const rates_func_options_base bopt(c);

    const unsigned n_branch_cat_sets(cm().branch_cat_set_size(c));
    for(unsigned bcs(0);bcs<n_branch_cat_sets;++bcs){

      std::string branch_set_label;
      cm().branch_cat_set_label(c,bcs,branch_set_label);

      const unsigned mmc(cm().typed_cat_no_from_cat_no_y_branch_cat_set(c,bcs,CAT_PARAM_TYPE::MUT_MODEL));
      const std::string& mmc_label(cm().typed_cat_label(mmc,CAT_PARAM_TYPE::MUT_MODEL));

      os << "\n";
      print_branch_cat_mut_model_header(bcs,branch_set_label,c,cat_label,mmc,mmc_label,os);

      if(rate_model_nuc(mmc) == RATE_MODEL_NUC::GY94){
        os << "GY94 nucleotide rates cannot be expressed in nucleotide tables\n"
           << "because of dependency on the codon distribution.\n\n";
        continue;
      }

      const CONTEXT_GROUP::index_t cg(context_model_group(context_model_nuc(mmc)));

      if(cg == CONTEXT_GROUP::TWOWAY){
        for(unsigned n0(0);n0<NUC::SIZE;++n0){
          for(unsigned n2(0);n2<NUC::SIZE;++n2){
            report_rate_set(static_cast<NUC::index_t>(n0),
                            static_cast<NUC::index_t>(n2),c,bcs,os);
          }
        }
        for(unsigned n0(0);n0<NUC::SIZE;++n0){
          NUC::index_t n2 = NUC::N;
          report_rate_set(static_cast<NUC::index_t>(n0),n2,c,bcs,os);
        }
        for(unsigned n2(0);n2<NUC::SIZE;++n2){
          NUC::index_t n0 = NUC::N;
          report_rate_set(n0,static_cast<NUC::index_t>(n2),c,bcs,os);
        }

      } else if(cg == CONTEXT_GROUP::ONLY_PRE){
        for(unsigned n0(0);n0<NUC::SIZE;++n0){
          NUC::index_t n2 = NUC::N;
          report_rate_set(static_cast<NUC::index_t>(n0),n2,c,bcs,os);
        }

      } else if(cg == CONTEXT_GROUP::ONLY_POST){
        for(unsigned n2(0);n2<NUC::SIZE;++n2){
          NUC::index_t n0 = NUC::N;
          report_rate_set(n0,static_cast<NUC::index_t>(n2),c,bcs,os);
        }
      }

      report_rate_set(NUC::N,NUC::N,c,bcs,os);


      if(model_type()==RATE_GTOR_MODEL::C4POST && n_cats==1 && cm().branch_cat_set_size(0)==1) {

        // get CpG/Non-CpG transition rate: start out by getting the state neutral flux matrix
        /// \todo (last minute hack!! -- fix me later)
        //
        const unsigned n_states(state_size());
        const unsigned N2(n_states*n_states);
        simple_array<smlfloat> flux(N2,0.);

        const prob_t* bgp(bg_pdistro_cat(0));

        rates(flux.ptr(),rates_func_options(rates_func_options_base(0,true,true),0));
        rates_to_flux_inplace(n_states,bgp,flux.ptr());

        smlfloat cpgti_flux(0.);
        smlfloat cpgti_p(0.);
        smlfloat noncpgti_flux(0.);
        smlfloat noncpgti_p(0.);
        smlfloat noncpgti_cg_flux(0.);
        smlfloat noncpgti_cg_p(0.);

        const SITE_MODEL::index_t sm(site_model());
        const unsigned base_size(SITE_MODEL::base_size(sm));

        NUC::index_t s1n[SITE_MODEL::MAX_BASE_SIZE];
        NUC::index_t s2n[SITE_MODEL::MAX_BASE_SIZE];

        for(unsigned s1(0);s1<n_states;++s1){
          SITE_MODEL::decode_nuc(sm,s1,s1n);

          for(unsigned s2(0);s2<n_states;++s2){
            SITE_MODEL::decode_nuc(sm,s2,s2n);

            // 1. is single nuc diff?
            unsigned nsubs(0);
            unsigned mut_p(0);
            for(unsigned p(0);p<base_size;++p){
              if( s1n[p] != s2n[p] ) {
                mut_p=p;
                nsubs++;
              }
            }

            if(nsubs != 1) continue;

            // 2. is it a transition?:
            if(! NUC::is_transition(s1n[mut_p],s2n[mut_p])) continue;

            // 3. is it a CpG,non-Cpg,or indeterminate?

            // don't use position 3 for anything else but the NGN->A
            if       (mut_p==3){
              if(s1n[mut_p] == NUC::G){
                if(s1n[mut_p-1] == NUC::C){
                  cpgti_flux += flux[s2+s1*n_states];
                  cpgti_p += bgp[s1];
                } else {
                  noncpgti_flux += flux[s2+s1*n_states];
                  noncpgti_p += bgp[s1];
                  noncpgti_cg_flux += flux[s2+s1*n_states];
                  noncpgti_cg_p += bgp[s1];
                }
              } else {
                //skip
              }

              // and dont read any G from pos 0
            } else if(mut_p==0 && s1n[mut_p] == NUC::G){
              //skip
            } else{
              if((s1n[mut_p] == NUC::C && s1n[mut_p+1] == NUC::G) ||
                 (s1n[mut_p] == NUC::G && s1n[mut_p-1] == NUC::C)){
                cpgti_flux += flux[s2+s1*n_states];
                cpgti_p += bgp[s1];
              } else {
                noncpgti_flux += flux[s2+s1*n_states];
                noncpgti_p += bgp[s1];
                if(s1n[mut_p] == NUC::C || s1n[mut_p] == NUC::G){
                  noncpgti_cg_flux += flux[s2+s1*n_states];
                noncpgti_cg_p += bgp[s1];
                }
              }
            }
          }
        }

        const smlfloat summ_cpgti(cpgti_flux/cpgti_p);
        const smlfloat summ_noncpgti(noncpgti_flux/noncpgti_p);
        const smlfloat summ_noncpg_cg_ti(noncpgti_cg_flux/noncpgti_cg_p);

        os << "SUMM_CPGTI: " << summ_cpgti << "\n";
        os << "SUMM_NONCPGTI: " << summ_noncpgti << "\n";
        os << "SUMM_CPGTI_RATIO: " << summ_cpgti/summ_noncpgti << "\n";
        os << "SUMM_NONCPGCGTI: " << summ_noncpg_cg_ti << "\n";
        os << "SUMM_CPGCGTI_RATIO: " << summ_cpgti/summ_noncpg_cg_ti << "\n";
        os << "\n";
        //       os << "TI: " << noncpgti_flux/noncpgti_p << "\n";
        //       os << "TV: " << noncpgti_flux/noncpgti_p << "\n";

      }
    }
  }
}



// i686 gcc 4.1 has lots of (harmless) link warnings when this is left inline:
//
void
rate_gtor_nuc_base::
bg_pdistro_update_internal(){
  base_t::bg_pdistro_update_internal();
  bg_pdistro_nuc_update();
  bg_pdistro_nuc_conditioned_update();
  _data.is_bg_valid=true;
}



void
rate_gtor_nuc_base::
set_mutation_model_slope(){
  static const smlfloat dull(0.9);

  const unsigned nmmc(cm().typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL));
  const smlfloat increment((1.-dull)*2./(nmmc+1));

  for(unsigned i(0);i<nmmc;++i){
    const unsigned st(paramset_mutation_model_cat_start(i));
    const unsigned s(paramset_mutation_model_cat_size(i));

    const smlfloat val(dull+(i+1)*increment);
    for(unsigned j(st);j<(st+s);++j) if(_is_train_param[j]) _param[j] *= val;
  }
}


CONTEXT_MODEL_NUC::index_t
rate_gtor_nuc_base::
context_model_nuc(const rates_func_options& opt) const {
  const unsigned smmc(cm().typed_cat_no_from_cat_no_y_branch_cat_set(opt.cat,opt.branch_cat_set,CAT_PARAM_TYPE::MUT_MODEL));

  return context_model_nuc(smmc);
}

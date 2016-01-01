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

// $Id: time_gtor.cc 1183 2008-03-27 02:12:11Z ctsa $

/// \file

#include "bi_tree.h"
#include "branch_time_util.h"
#include "cat_manager.h"
#include "subs_ml_print_util.h"
#include "substk_exception.h"
#include "time_gtor.h"
#include "util/general/io_util.h"
#include "util/general/log.h"
#include "util/math/random_util.h"
#include "util/math/test_float.h"

#include <iomanip>
#include <iostream>
#include <sstream>


void
time_gtor::
set_pseudo_param(const unsigned param_no,
                 const smlfloat val){
  check_pseudoparam_num(param_no);

  // pseudoparam n: scale all trainable times of cat n as a group:

  const smlfloat safe_val(std::fabs(val));
  const unsigned bs(tree().branch_size());
  for(unsigned i(0);i<bs;++i) {
    const unsigned ti(time_cat_index(i,param_no));
    if(_is_train_param[ti])  _param[ti] *= safe_val;
  }
}



smlfloat
time_gtor::
trainable_param_transform(const unsigned i) const {
  // time is presented to the minimizer in a more robust format,
  // meant to prevent "popups" to very high time values

  // adjust down from zero thresh slightly so that we can be
  // sure that x = exp(log(x)) is still <= zero thresh
  static const smlfloat fudge_down(0.99);
  smlfloat time_tmp(std::max(_param[i],ZERO_BRANCH_TIME_THRESHOLD*fudge_down));
  if(time_tmp>1.) time_tmp *= time_tmp;
  return std::log(time_tmp);
}



void
time_gtor::
trainable_param_untransform(const unsigned i,
                            const smlfloat p){
  smlfloat time_tmp(std::exp(p));
  if(time_tmp>1.) time_tmp = std::sqrt(time_tmp);
  _param[i] = time_tmp;
}



void
time_gtor::
post_write_param_state() {
  _time_penalty=0.;
  const unsigned ps(param_size());
  for(unsigned i(0);i<ps;++i){
    if(_is_train_param[i]) adjust_time(_param[i]);
  }

  check_param();
}



void
time_gtor::
adjust_time(smlfloat& time) {
  static const smlfloat MAX_BRANCH_TIME(1.5);
  if(time > MAX_BRANCH_TIME) {
    _time_penalty += (time - MAX_BRANCH_TIME);
    log_os  << "time penalty triggered, converted from: " << time << " to " << MAX_BRANCH_TIME << "\n";
    time = MAX_BRANCH_TIME;
  } else if(time<=ZERO_BRANCH_TIME_THRESHOLD) {
    time = 0.;
  }
}



const char* const section_id = "time_gtor";
const char* const end_label = "END";

const iliner il(section_id);



void
time_gtor::
store_state(std::ostream& os) const {

  const unsigned ntc(time_cat_size());
  const unsigned bs(tree().branch_size());

  os << std::setprecision(SUBS_ML_PRINT_MODEL_PRECISION);
  for(unsigned c(0);c<ntc;++c){
    for(unsigned b(0);b<bs;++b){
      const unsigned ci(time_cat_index(b,c));
      os << section_id << " time_cat_"<< c << "_branch_" << b << " " << tree().branch_node(b)->label() << " " <<
        _is_train_param[ci] << " " <<
        _param[ci] << "\n";
    }
  }

  os << section_id << " " << end_label << "\n";
}



time_gtor::
time_gtor(std::istream& is,
          const time_gtor_sml_share& smls)
    : base_t(), _smls(smls), _time_penalty(0) {

  std::string dummy;

  setup_param();

  const unsigned bs(tree().branch_size());
  const unsigned ntc(time_cat_size());

  for(unsigned c(0);c<ntc;++c){
    for(unsigned i(0);i<bs;++i){
      il.advance(is); // time_cat_i
      is >> dummy;
      const unsigned bi(tree().node(dummy)->branch_id());
      const unsigned ci(time_cat_index(bi,c));
      bool is_t;
      is >> is_t >> _param[ci];
      _is_train_param[ci] = is_t;
    }
  }

  check_param();

  il.advance(is,end_label);
}




void
time_gtor::
reset_param(const PARAM_INIT_TYPE::index_t pinit){

  if(pinit == PARAM_INIT_TYPE::RANDOM){
    static const smlfloat LN_TEN(std::log(10.));
    static const smlfloat time_min(0.005);
    static const smlfloat time_max(0.25);
    static const smlfloat log10_time_min(std::log(time_min)/LN_TEN);
    static const smlfloat log10_time_max(std::log(time_max)/LN_TEN);
    static const smlfloat timerange(log10_time_max-log10_time_min);

    const unsigned ps(param_size());
    for(unsigned i(0);i<ps;++i){
      if(_is_train_param[i]){
        const smlfloat rand_log10(random_uniform()*timerange+log10_time_min);
        _param[i] = std::exp(rand_log10*LN_TEN);
      }
    }

    tree_nonconst().clear_branch_lengths();
  } else if(pinit == PARAM_INIT_TYPE::START){
    static const smlfloat start_time(0.02);

    /// \todo give each time category a little init time separation,
    /// so that they can be pulled apart from the init state by a
    /// gradient minmizer:
    ///
    const unsigned ps(param_size());
    for(unsigned i(0);i<ps;++i){
      if(_is_train_param[i]){
        _param[i] = start_time;
      }
    }

    // don't setup initial category parameters on a slope if there is
    // only one per adset:
    if(cm().assigned_data_set_size() != 1){ die("finish seting up time init for multi-adset"); }
    set_time_cat_slope();

    // pull init times from tree if:
    // 1) they already existed in the newick string
    // 2) there is only one time category
    //
    bool is_warn_clear(true);
    if(time_cat_size() == 1){
      for(unsigned i(0);i<ps;++i){
        if(_is_train_param[i]){
          if(tree().branch_node(i)->is_length_set()){
            _param[i] = tree().branch_node(i)->length();
          }
        }
      }
      is_warn_clear=false;
    }

    tree_nonconst().clear_branch_lengths(is_warn_clear);

  } else { die("unknown param_init_type"); }
}



void
time_gtor::
check_param() const {
  const unsigned ps(param_size());
  for(unsigned i(0);i<ps;++i) {
    if(is_float_invalid(_param[i])){
      std::ostringstream oss;
      oss << "time_gtor: Invalid model parameter: " << i << " " << _param[i] << "\n";
      log_os << oss.str();
      throw substk_exception(oss.str().c_str());
    }
  }

  for(unsigned i(0);i<ps;++i){
    if(_param[i]<0.){
      std::ostringstream oss;
      oss << "time_gtor: Negative branch time parameter: " << i << " " << _param[i] << "\n";
      log_os << oss.str();
      throw substk_exception(oss.str().c_str());
    }
  }
}



void
time_gtor::
dump(std::ostream& os) const {
  const unsigned ntc(time_cat_size());
  const unsigned bs(tree().branch_size());
  for(unsigned c(0);c<ntc;++c){
    for(unsigned i(0);i<bs;++i){
      const unsigned ci(time_cat_index(i,c));
      os << tree().branch_node(i)->label() << "\t" << _param[ci] << "\n";
    }
  }
}



void
time_gtor::
set_reversible_time_param(const unsigned time_cat){

  const bi_tree& t(tree());

  if(t.root()->is_leaf()) return;
  const bi_tree_node* c1(t.root()->child1());
  const bi_tree_node* c2(t.root()->child2());

  // check if either branch is already locked:
  const unsigned bi1(c1->branch_id());
  const unsigned pi1(time_cat_index(bi1,time_cat));

  const unsigned bi2(c2->branch_id());
  const unsigned pi2(time_cat_index(bi2,time_cat));
  if( _is_train_param[pi1]==false || _is_train_param[pi2]==false) return;

  // lock child1, unless it's a leaf node:
  unsigned pix(pi1);
  if(c1->is_leaf()) pix=pi2;

  _param[pix]=0.;
  _is_train_param[pix]=false;
}



void
time_gtor::
set_time_cat_slope(){
  static const smlfloat dull(0.9);

  const unsigned ntc(time_cat_size());
  const smlfloat increment((1.-dull)*2./(ntc+1));
  const unsigned bs(tree().branch_size());
  for(unsigned c(0);c<ntc;++c){
    const smlfloat val(dull+(c+1)*increment);
    for(unsigned i(0);i<bs;++i){
      const unsigned ci(time_cat_index(i,c));
      if(_is_train_param[ci]) _param[ci] *= val;
    }
  }
}



bool
time_gtor::
is_zero_time(const unsigned cat) const {
  const unsigned bs(tree().branch_size());
  smlfloat branch_sum(0.);
  for(unsigned i(0);i<bs;++i){
    branch_sum += branch_time(i,cat);
  }
  return (branch_sum <= 0.);
}



unsigned
time_gtor::
time_cat_index(const unsigned i,
               const unsigned time_cat) const {
  const unsigned n_branches(tree().branch_size());
  assert(i<n_branches);
  return time_cat*n_branches+i;
}



void
time_gtor::
setup_param(){
  resize(tree().branch_size()*time_cat_size());
  const unsigned ps(param_size());
  for(unsigned i(0);i<ps;++i){
    _is_nonnegative[i] = false;
  }
}



unsigned
time_gtor::
time_cat_size() const {
  return cm().typed_cat_size(CAT_PARAM_TYPE::TIME);
}



void
time_gtor::
time_cat_pdistro(smlfloat* p) const {
  cm().typed_cat_pdistro(p,CAT_PARAM_TYPE::TIME);
}



unsigned
time_gtor::
time_cat_no(const unsigned cat) const {
  return cm().typed_cat_no_from_cat_no(cat,CAT_PARAM_TYPE::TIME);
}



smlfloat
time_gtor::
cat_averaged_branch_time(const unsigned i) const {
  const unsigned n_time_cats(time_cat_size());

  prob_t val(0.);
  for(unsigned tc(0);tc<n_time_cats;++tc){
    val += time_cat_branch_time(i,tc);
  }
  return val;
}



smlfloat
time_gtor::
cat_expected_branch_time(const unsigned i) const {
  const unsigned n_time_cats(time_cat_size());
  simple_array<prob_t> tcp(n_time_cats);
  time_cat_pdistro(tcp.ptr());
  prob_t val(0.);
  for(unsigned tc(0);tc<n_time_cats;++tc){
    val += time_cat_branch_time(i,tc)*tcp[tc];
    }
  return val;
}

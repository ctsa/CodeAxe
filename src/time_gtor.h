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

// $Id: time_gtor.h 1077 2008-01-05 05:00:19Z ctsa $

/// \file

#ifndef __TIME_GTOR_H
#define __TIME_GTOR_H

#include "param_init_type.h"
#include "param_object.h"
#include "simple_util.h"
#include "subs_ml_model_comex.h"
#include "subs_ml_types.h"
#include "util/general/uncopyable.h"
#include "util/math/indy_random_var.h"

#include <string>

struct bi_tree;
struct cat_manager;



struct time_gtor_options {

  time_gtor_options() : is_reversible(false) {}

  bool is_reversible;
};



struct time_gtor : public param_object, private uncopyable {

  typedef param_object base_t;
  typedef time_gtor self_t;

  time_gtor(std::istream& is,
            const time_gtor_sml_share& smls);

  time_gtor(const time_gtor_options& tgo,
            const time_gtor_sml_share& smls)
    : base_t(),
      _smls(smls),
      _time_penalty(0.) {
    setup_param();
    if(tgo.is_reversible) set_reversible_time_param();
  }

  time_gtor(const self_t& s,
            const time_gtor_sml_share& smls)
    : base_t(s),
      _smls(smls),
      _time_penalty(s._time_penalty) {}

  smlfloat
  branch_time(const unsigned i,
              const unsigned cat) const {
    return _param[cat_index(i,cat)];
  }

  smlfloat
  time_cat_branch_time(const unsigned i,
                       const unsigned time_cat) const {
    return _param[time_cat_index(i,time_cat)];
  }

  smlfloat
  cat_averaged_branch_time(const unsigned i) const;

  smlfloat
  cat_expected_branch_time(const unsigned i) const;

  /// \todo this isn't cool... it relies on the index ordering, think
  /// about gitten' rid of it
  ///
  std::vector<smlfloat>::const_iterator branch_times(const unsigned cat) const {
    return _param.begin()+cat_index(0,cat);
  }

  const irv_t<smlfloat> branch_time_variance(const unsigned i,
                                             const unsigned cat) const {
    const unsigned index(cat_index(i,cat));
    return irv_t<smlfloat>(_param[index],param_variance(index));
  }

  const irv_t<smlfloat> time_cat_branch_time_variance(const unsigned i,
                                                      const unsigned time_cat) const {
    const unsigned index(time_cat_index(i,time_cat));
    return irv_t<smlfloat>(_param[index],param_variance(index));
  }

  void store_state(std::ostream& os) const;

  /// \brief quick debugging output of model parameters
  void dump(std::ostream& os) const;

#if 0
  /// \brief report information processed from model parameters
  void report(std::ostream& os) const;

  /// \brief time subsection of report
  void report_time(std::ostream& os) const;
#endif

  const smlfloat& param_penalty() const { return _time_penalty; }

  void reset_param(const PARAM_INIT_TYPE::index_t pinit);

  /// client-code can decide what to do about this case:
  ///
  bool is_zero_time(const unsigned cat) const;

  unsigned time_cat_size() const;

  void time_cat_pdistro(smlfloat* p) const;

#if 0
  void
  time_cat_time_scale(const unsigned time_cat,
                      const smlfloat scale) {
    const unsigned bs(branch_size());
    for(unsigned i(0);i<bs;++i) _param[time_cat_index(i,time_cat)] *= scale;
  }
#endif

private:
  const bi_tree& tree() const { return _smls.tree(); }
  bi_tree& tree_nonconst() { return _smls.tree_nonconst(); }

  /// adjusts all time parameters for a reversible model
  ///
  void set_reversible_time_param(){
    const unsigned n_time_cats(time_cat_size());
    for(unsigned tc(0);tc<n_time_cats;++tc){
      set_reversible_time_param(tc);
    }
  }

  /// adjusts time category parameters for a reversible model
  ///
  void set_reversible_time_param(const unsigned time_cat);

  unsigned time_cat_no(const unsigned cat) const;

  unsigned
  cat_index(const unsigned i,
            const unsigned cat) const {
    return time_cat_index(i,time_cat_no(cat));
  }

  unsigned
  time_cat_index(const unsigned i,
                 const unsigned time_cat) const;

  const cat_manager&
  cm() const {
    return _smls.get_cat_manager();
  }

  virtual
  unsigned pseudo_param_size() const { return time_cat_size(); } // pseudoparam to scale all times together

  virtual
  void
  set_pseudo_param(const unsigned param_no,
                   const smlfloat val);

  virtual
  smlfloat
  trainable_param_transform(const unsigned i) const;

  virtual
  void
  trainable_param_untransform(const unsigned i,
                              const smlfloat p);

  virtual
  void
  post_write_param_state();

  void check_param() const;

  void adjust_time(smlfloat& time);

  void check_pseudoparam_num(unsigned param_no) const {
    if(param_no>=pseudo_param_size()){
      die("Bad pseudoparam number");
    }
  }

  void setup_param();

  void set_time_cat_slope();

  //////
  time_gtor_sml_share _smls;
  smlfloat _time_penalty;
};

#endif

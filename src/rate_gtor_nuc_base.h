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

// $Id: rate_gtor_nuc_base.h 1133 2008-01-29 22:42:35Z ctsa $

/// \file

#ifndef __RATE_GTOR_NUC_BASE_H
#define __RATE_GTOR_NUC_BASE_H

#include "cat_type.h"
#include "rate_gtor.h"
#include "rate_gtor_nuc_options.h"
#include "simple_util.h"
#include "util/bio/bioseq_util.h"
#include "util/math/indy_random_var.h"

#include <cassert>
#include <cmath>
#include <iosfwd>


/// \brief base rate generator for all models using nuc mutation
/// (coding or not)
///
struct rate_gtor_nuc_base : public rate_gtor {
  typedef rate_gtor base_t;
  typedef rate_gtor_nuc_base self_t;

  rate_gtor_nuc_base(const rate_gtor_options& ropt,
                     const rate_gtor_nuc_options& nopt,
                     const rate_gtor_sml_share& smls)
    : base_t(ropt,smls), _data() {
    _data.nopt = nopt;
    init_class_storage();
    init_param_locks(nopt.lopt);
    if(! nopt.is_free_edge_correction) lock_full_edge_correction();
  }

  rate_gtor_nuc_base(const self_t& s,
                     const rate_gtor_sml_share& smls)
    :  base_t(s,smls), _data(s._data) {}

  virtual
  void rates(smlfloat* rates,
             const rates_func_options& opt) const;

  virtual
  void store_state(std::ostream& os) const;

  virtual
  void report(std::ostream& os) const;

  smlfloat
  cat_nuc_context_mut_rate(const NUC::index_t from0,
                           const NUC::index_t from1,
                           const NUC::index_t from2,
                           const NUC::index_t to,
                           const unsigned cat_no,
                           const unsigned branch_cat_set,
                           const prob_t* nuc_bg0,
                           const prob_t* nuc_bg2,
                           const bool is_coding_model = false) const;

  const prob_t* bg_pdistro_nuc_cat(const unsigned cat_no) const {
    bg_check();
    return  _data.bg_pdistro_nuc_bg_cat[bg_cat_no_from_cat_no(cat_no)];
  }

  const prob_t*
  bg_pdistro_nuc_conditioned_on_5p_cat(const unsigned cat_no,
                                       const NUC::index_t dep_n = static_cast<NUC::index_t>(0)) const{
    bg_check();
    return _data.bg_pdistro_nuc_conditioned_on_5p_bg_cat[bg_cat_no_from_cat_no(cat_no)]+dep_n*NUC::SIZE;
  }

  const prob_t*
  bg_pdistro_nuc_conditioned_on_3p_cat(const unsigned cat_no,
                                       const NUC::index_t dep_n = static_cast<NUC::index_t>(0)) const{
    bg_check();
    return _data.bg_pdistro_nuc_conditioned_on_3p_bg_cat[bg_cat_no_from_cat_no(cat_no)]+dep_n*NUC::SIZE;
  }

  CONTEXT_MODEL_NUC::index_t
  context_model_nuc(const unsigned i) const { return _data.nopt.context_model(i); }

  CONTEXT_MODEL_NUC::index_t
  context_model_nuc(const rates_func_options& opt) const;

  bool
  is_use_edge_correction() const { return _data.nopt.is_use_edge_correction; }

  smlfloat
  edge_strength_pre() const { return edge_strength(0); }

  smlfloat
  edge_strength_post() const { return edge_strength(1); }
#if 1
  smlfloat
  edge_strength_pre2() const { return edge_strength(0); }

  smlfloat
  edge_strength_post2() const { return edge_strength(1); }
#endif

  unsigned
  param_size_mut() const {
    return paramset_size(PARAMSET_NUC);
  }

  unsigned
  param_start_mut() const {
    return paramset_start(PARAMSET_NUC);
  }

  ///////////////////////////////////////////////////////////
protected:

  virtual
  void reset_param_internal(const PARAM_INIT_TYPE::index_t pinit);

  virtual void bg_pdistro_update_internal();

  virtual void load_state_internal(std::istream& is);

  void stat_pdistro_nuc(prob_t* nuc_sd,
                        const rates_func_options& opt) const;

  smlfloat
  edge_strength(const unsigned i) const {
    if(! is_use_edge_correction()) return 0.;
    assert(i<paramset_size(PARAMSET_EDGE_STRENGTH));
    return _param[paramset_start(PARAMSET_EDGE_STRENGTH)+i];
  }

  unsigned
  param_end() const { return paramset_start(PARAMSET_END); }
#if 0
  // cat sizes:
  unsigned
  cat_size_site_rate() const { return cm().typed_cat_size(CAT_PARAM_TYPE::MUT_RATE,CAT_MIX_TYPE::SITE); }

  unsigned
  cat_size_mutation_model() const { return cm().typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL); }

  unsigned
  cat_size_group_rate() const { return  cm().typed_cat_size(CAT_PARAM_TYPE::MUT_RATE,CAT_MIX_TYPE::GROUP);}
#endif
  virtual void fix_input_param_internal();

  virtual void param_update_from_new_bg_pdistro_internal();

  virtual
  void state_pdistro_to_nuc_pdistro(const prob_t* state_pdistro,
                                    prob_t* nuc_pdistro) const = 0;

  virtual
  smlfloat
  trainable_param_transform(const unsigned i) const {
    if(i<base_t::param_end()) {
      return base_t::trainable_param_transform(i);
    }
    const unsigned st(paramset_start(PARAMSET_EDGE_STRENGTH));
    const unsigned s(paramset_size(PARAMSET_EDGE_STRENGTH));
    if(i>=st && i<(st+s)){
      return std::acos((_param[i]*2.)-1.);
    } else {
      return _param[i];
    }
  }

  virtual
  void
  trainable_param_untransform(const unsigned i,
                              const smlfloat p) {
    if(i<base_t::param_end()) {
      return base_t::trainable_param_untransform(i,p);
    }
    const unsigned st(paramset_start(PARAMSET_EDGE_STRENGTH));
    const unsigned s(paramset_size(PARAMSET_EDGE_STRENGTH));
    if(i>=st && i<(st+s)){
      _param[i] = (std::cos(p)+1.)/2.;
    } else {
      _param[i] = p;
    }
  }

  // parameter arrangement:
  enum paramset_t {
    PARAMSET_SITE_RATE_CATS,
    PARAMSET_GROUP_RATE_CATS,
    PARAMSET_NUC,
    PARAMSET_EDGE_STRENGTH,
    PARAMSET_END
  };

  virtual
  unsigned
  paramset_size(const int p) const;

  unsigned
  paramset_start(const int p) const {
    unsigned s(base_t::param_end());
    for(int i(0);i<p;++i) s+=paramset_size(i);
    return s;
  }

  void
  set_paramset_dependent_block(const int p){
    const unsigned st(paramset_start(p));
    const unsigned s(paramset_size(p));
    if(s){
      _is_free_block_start[st] = true;
      _is_free_block_stop[st+s-1] = true;
    }
  }

  void
  set_paramset_val(const int p,
                   const smlfloat val){
    const unsigned st(paramset_start(p));
    const unsigned s(paramset_size(p));
    for(unsigned i(st);i<(st+s);++i)  _param[i] = val;
  }

  /// set paramset values in even increments from 2. to 0.  this is
  /// intended to minimally pull initial values apart so that each has
  /// a distinct gradient during initial minimization
  ///
  void
  set_paramset_slope(const int p){
    static const smlfloat dull(0.9);

    const unsigned st(paramset_start(p));
    const unsigned s(paramset_size(p));
    const smlfloat increment((1.-dull)*2./(s+1));
    for(unsigned i(0);i<(s);++i) {
      if(_is_train_param[i+st]) _param[i+st] = dull+(i+1)*increment;
    }
  }

  void
  set_paramset_train(const int p,
                     bool b){
    const unsigned st(paramset_start(p));
    const unsigned s(paramset_size(p));
    for(unsigned i(st);i<(st+s);++i) _is_train_param[i]=b;
  }

  irv_t<smlfloat>
  get_paramset_member(const int p,
                      const unsigned i) const {
    assert(i<paramset_size(p));
    const unsigned pi(paramset_start(p)+i);
    return irv_t<smlfloat>(_param[pi],param_variance(pi));
  }

  RATE_MODEL_NUC::index_t
  rate_model_nuc(const unsigned i) const { return _data.nopt.rate_model(i); }

  void
  fix_input_param_cat_branch_paramset(const int paramset,
                                      const CAT_PARAM_TYPE::index_t pt,
                                      const CAT_MIX_TYPE::index_t mt,
                                      const unsigned n_sets,
                                      const std::vector<unsigned>::const_iterator& set_from_cat);

  void
  lock_1param_norm_sets(const int paramset,
                        const CAT_PARAM_TYPE::index_t pt,
                        const CAT_MIX_TYPE::index_t mt,
                        const unsigned n_sets,
                        const std::vector<unsigned>::const_iterator set_from_cat);

  ///////////////////////////////////////////////////////////////////
private:
  void init_class_storage();

  void init_param_locks(const rate_gtor_nuc_lock_options& lopt);

  void
  paramset_label(const paramset_t p,
                 unsigned n,
                 std::ostream& os) const;

  void
  scale_nuc_mut_rate_trainable_param(const unsigned mut_model_cat,
                                     smlfloat x){
    const unsigned st(paramset_mutation_model_cat_start(mut_model_cat));
    const unsigned s(paramset_mutation_model_cat_size(mut_model_cat));
    if(context_model_nuc(mut_model_cat) == CONTEXT_MODEL_NUC::FACTORED_TRIPLET) x = std::sqrt(x);
    for(unsigned i(st);i<(st+s);++i) if(_is_train_param[i]) _param[i] *= x;
  }

  struct trainable_mode {
    enum index_t { NONE, FULL, HALF_JOINT };
  private:
    trainable_mode() {}
  };


  /// Internally used to scale model file parameters.
  ///
  /// Reports the time scale of the rate matrix for one mutation model
  /// category using the context-independent bg nucleotide
  /// distribution associated with cat_no, averaging over all edge
  /// contexts
  ///
  /// trainable mode facilitates normalization when some rate
  /// parameters are locked
  ///
  smlfloat rate_scale_cat_no_y_branch_id(const unsigned cat_no,
                                         const unsigned branch_id,
                                         const trainable_mode::index_t tm = trainable_mode::NONE) const;

  enum joint_param_t { JOINT1, JOINT2 };

  void
  get_nuc_param_index_name(const unsigned mut_model_param_no,
                           const unsigned param_index,
                           std::ostream& os) const;

  unsigned
  nuc_context_mut_rate_param_index(const NUC::index_t from0,
                                   NUC::index_t from1,
                                   const NUC::index_t from2,
                                   NUC::index_t to,
                                   const unsigned mut_model_cat,
                                   const joint_param_t jp = JOINT1) const;

  irv_t<smlfloat>
  nuc_context_mut_rate_param(const NUC::index_t from0,
                             const NUC::index_t from1,
                             const NUC::index_t from2,
                             const NUC::index_t to,
                             const unsigned mut_model_cat,
                             const trainable_mode::index_t tm) const;

  irv_t<smlfloat>
  nuc_context_mut_rate(const NUC::index_t from0,
                       const NUC::index_t from1,
                       const NUC::index_t from2,
                       const NUC::index_t to,
                       const unsigned cat_no,
                       const unsigned branch_cat_set,
                       const prob_t* nuc_bg0,
                       const prob_t* nuc_bg2,
                       const trainable_mode::index_t tm = trainable_mode::NONE) const;


  // init to a fixed transition/transversion ratio
  //
  void
  set_ts_tv_ratio(const unsigned mut_model_cat,
                  const smlfloat ratio);


  // helper func to set_ts_tv_ratio:
  void
  set_ts_tv_ratio_nhood(const unsigned mut_model_cat,
                        const smlfloat ratio,
                        const NUC::index_t n0,
                        const NUC::index_t n2);

  // helper func to report
  void
  report_rate_set(NUC::index_t n0,
                  NUC::index_t n2,
                  const unsigned cat_no,
                  const unsigned branch_cat_set,
                  std::ostream& os) const;

  void
  report_rate_cats(const CAT_MIX_TYPE::index_t mt,
                   const int paramset,
                   const unsigned branch_id,
                   std::ostream& os) const;

  static
  unsigned
  get_rate_context_model_param_size(const RATE_MODEL_NUC::index_t rm,
                                    const CONTEXT_MODEL_NUC::index_t cm);

  unsigned
  paramset_mutation_model_cat_size(const unsigned cat_no) const {
    return get_rate_context_model_param_size(rate_model_nuc(cat_no),
                                             context_model_nuc(cat_no));
  }

  unsigned
  paramset_mutation_model_cat_start(const unsigned cat_no) const {
    unsigned val(paramset_start(PARAMSET_NUC));

    for(unsigned i(0);i<cat_no;++i){
      val += paramset_mutation_model_cat_size(i);
    }

    return val;
  }

  void set_mutation_model_slope();

  prob_t* bg_pdistro_nuc_cat(const unsigned cat_no) {
    return bg_pdistro_nuc_bg_cat(bg_cat_no_from_cat_no(cat_no));
  }

  prob_t* bg_pdistro_nuc_bg_cat(const unsigned bg_cat_no) {
    return _data.bg_pdistro_nuc_bg_cat[bg_cat_no];
  }

  void
  bg_pdistro_nuc_update(){
    const unsigned dcn(bg_cat_size());
    _data.bg_pdistro_nuc_bg_cat.init(dcn,NUC::SIZE);
    for(unsigned c(0);c<dcn;++c){
      state_pdistro_to_nuc_pdistro(bg_pdistro_bg_cat(c),bg_pdistro_nuc_bg_cat(c));
    }
  }

  void bg_pdistro_nuc_conditioned_update();

  void bg_check() const {
    if(!_data.is_bg_valid) die("rate_gtor_nuc_base: invalid pdistro bg state");
  }

  void normalize_mut_rate_params();

  void
  lock_full_edge_correction() {
    const unsigned st(paramset_start(PARAMSET_EDGE_STRENGTH));
    const unsigned s(paramset_size(PARAMSET_EDGE_STRENGTH));
    for(unsigned i(st);i<(st+s);++i){
      _param[i] = 1.;
      _is_train_param[i]=false;
    }
  }

  /// struct to simplify class copy constructor:
  struct auto_copy {
    auto_copy()
      : is_bg_valid(false) {}

    rate_gtor_nuc_options_base nopt;

    simple_init_matrix<prob_t> bg_pdistro_nuc_bg_cat;
    simple_init_matrix<prob_t> bg_pdistro_nuc_conditioned_on_5p_bg_cat;
    simple_init_matrix<prob_t> bg_pdistro_nuc_conditioned_on_3p_bg_cat;

    bool is_bg_valid;
  };

  ///////////////////////////data:
  auto_copy _data;
};


#endif

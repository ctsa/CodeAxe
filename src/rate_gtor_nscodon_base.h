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

// $Id: rate_gtor_nscodon_base.h 1146 2008-03-01 00:21:53Z ctsa $

/// \file

#ifndef __RATE_GTOR_NSCODON_BASE_H
#define __RATE_GTOR_NSCODON_BASE_H


#include "rate_gtor_nuc_base.h"
#include "rate_gtor_nscodon_options.h"

#include <iosfwd>
#include <utility>

namespace AA_PARAM {
  enum index_t {
    normal,
    zero,
    one };
}

/// \brief base rate generator for coding sequence state models
///
struct rate_gtor_nscodon_base : public rate_gtor_nuc_base {

  typedef rate_gtor_nuc_base base_t;
  typedef rate_gtor_nscodon_base self_t;

  rate_gtor_nscodon_base(const rate_gtor_options& ropt,
                         const rate_gtor_nuc_options& nopt,
                         const rate_gtor_nscodon_options& copt,
                         const rate_gtor_sml_share& smls)
    : base_t(ropt,nopt,smls), _data() {
    _data.copt = copt;
    init_class_storage();
    init_param_locks(copt.lopt);
  }

  rate_gtor_nscodon_base(const self_t& s,
                         const rate_gtor_sml_share& smls)
  : base_t(s,smls), _data(s._data) {}

  virtual
  void rates(smlfloat rates[],
             const rates_func_options& opt) const;

  virtual
  void rates_variance(irv_t<smlfloat> rates[],
                      const rates_func_options& opt) const;

  virtual
  void store_state(std::ostream& os) const;

  /// \brief prints nuc rates and aa selection parameters
  virtual
  void report(std::ostream& os) const;

  irv_t<smlfloat>
  cat_nsaa_param(NSAA::index_t from,
                 NSAA::index_t to,
                 const unsigned cat_no,
                 const unsigned branch_cat_set) const;

  smlfloat
  codon_bias(const NSCODON::index_t from,
             const NSCODON::index_t to) const;
#if 0
  bool
  is_codon_bias() const {
    return (_data.copt.codon_bias_model != CODON_BIAS_MODEL::NONE);
  }
#endif

  const prob_t* bg_pdistro_nuc_pos_cat(const unsigned cat_no,
                                       const unsigned pos) const {
    bg_check();
    return  _data.bg_pdistro_nuc_pos_bg_cat[bg_cat_no_from_cat_no(cat_no)][pos];
  }

  const prob_t*
  bg_pdistro_nuc_pos_conditioned_on_5p_cat(const unsigned cat_no,
                                           const unsigned pos,
                                           const NUC::index_t dep_n = static_cast<NUC::index_t>(0)) const{
    bg_check();
    return _data.bg_pdistro_nuc_pos_conditioned_on_5p_bg_cat[bg_cat_no_from_cat_no(cat_no)][pos]+dep_n*NUC::SIZE;
  }

  const prob_t*
  bg_pdistro_nuc_pos_conditioned_on_3p_cat(const unsigned cat_no,
                                           const unsigned pos,
                                           const NUC::index_t dep_n = static_cast<NUC::index_t>(0)) const{
    bg_check();
    return _data.bg_pdistro_nuc_pos_conditioned_on_3p_bg_cat[bg_cat_no_from_cat_no(cat_no)][pos]+dep_n*NUC::SIZE;
  }

  const prob_t* bg_pdistro_nscodon_cat(const unsigned cat_no) const {
    bg_check();
    return  _data.bg_pdistro_nscodon_bg_cat[bg_cat_no_from_cat_no(cat_no)];
  }

  virtual
  smlfloat param_penalty() const;

  // public versions of these methods to set ci subset (temporary!!)
  unsigned
  param_size_aa() const {
    return paramset_size(PARAMSET_AA);
  }

  unsigned
  param_start_aa() const {
    return paramset_start(PARAMSET_AA);
  }

  smlfloat
  cat_nuc_context_mut_rate(const NUC::index_t from0,
                           const NUC::index_t from1,
                           const NUC::index_t from2,
                           const NUC::index_t to,
                           const unsigned cat_no,
                           const unsigned branch_cat_set,
                           const prob_t* nuc_bg0,
                           const prob_t* nuc_bg2,
                           const NSCODON::index_t to_state) const;

  /////////////////////////////////////////////////////////////////////////////
protected:

  virtual
  void load_state_internal(std::istream& is);

  virtual
  void reset_param_internal(const PARAM_INIT_TYPE::index_t pinit);

  virtual void bg_pdistro_update_internal();

  unsigned
  param_end() const { return paramset_start(PARAMSET_END); }

  virtual
  unsigned pseudo_param_size() const;

  virtual
  void set_pseudo_param(const unsigned param_no,
                        const smlfloat val);

  std::pair<unsigned,AA_PARAM::index_t>
  nsaa_param_index(NSAA::index_t from,
                   NSAA::index_t to,
                   const unsigned select_matrix_cat) const;

  irv_t<smlfloat>
  nsaa_param(const NSAA::index_t from,
             const NSAA::index_t to,
             const unsigned select_matrix_cat) const;

  void state_pdistro_to_nscodon_pdistro(const prob_t* state_pdistro,
                                        prob_t* nscodon_pdistro) const;

  void state_pdistro_to_nsaa_pdistro(const prob_t* state_pdistro,
                                     prob_t* nsaa_pdistro) const;

  // parameter arrangement:
  enum paramset_t {
    PARAMSET_SITE_SELECT_CATS = base_t::PARAMSET_END,
    PARAMSET_GROUP_SELECT_CATS,
    PARAMSET_AA,
    PARAMSET_CODON_BIAS,
    PARAMSET_END
  };

  virtual
  unsigned paramset_size(const int p) const;

  /////////////////////////////////////////////////////////////////////////////
private:
  virtual
  void fix_input_param_internal();

  virtual
  void state_pdistro_to_nuc_pdistro(const prob_t* state_pdistro,
                                    prob_t* nuc_pdistro) const;

  void stat_pdistro_nscodon(prob_t* sd,
                            const rates_func_options& opt) const;

  void stat_pdistro_nsaa(prob_t* sd,
                         const rates_func_options& opt) const;

  void bg_pdistro_nscodon_update(){
    const unsigned bcn(bg_cat_size());
    _data.bg_pdistro_nscodon_bg_cat.init(bcn,NSCODON::SIZE);
    for(unsigned c(0);c<bcn;++c){
      state_pdistro_to_nscodon_pdistro(bg_pdistro_bg_cat(c),bg_pdistro_nscodon_bg_cat(c));
    }
  }

  void bg_pdistro_nuc_pos_update();

  void bg_pdistro_nuc_pos_conditioned_update();

  void init_class_storage();

  void init_param_locks(const rate_gtor_nscodon_lock_options& lopt);

  // dumps a string identifying each paramset to os
  void
  paramset_label(const paramset_t p,
                 unsigned n,
                 std::ostream& os) const;

  void
  get_aa_param_index_name(const unsigned matrix_cat_no,
                          const unsigned matrix_param_no,
                          std::ostream& os) const;

  /// other:
  SELECT_MODEL::index_t select_model(const unsigned i) const { return _data.copt.select_model(i); }

  void
  nsaa_selection_matrix(irv_t<smlfloat> nsaa_selection[NSAA::SIZE*NSAA::SIZE],
                        const unsigned cat_no,
                        const unsigned branch_cat_set) const;

  smlfloat
  codon_bias_factor(const NSCODON::index_t to) const;

  /// \todo shouldn't both of these automatically have an expect of 1 b/c of the
  /// fix_input_param() call? both of these functions should go...
  ///
  smlfloat
  site_select_cats_expect() const {
#if 0
    const unsigned st = paramset_start(PARAMSET_SITE_SELECT_CATS);
    const unsigned s = paramset_size(PARAMSET_SITE_SELECT_CATS);

    smlfloat sum(0.);
    for(unsigned i(0);i<(s);++i){ sum += _param[st+i]*get_paramset_pdistro_member(PARAMSET_SITE_SELECT_CATS_PROB,i); }
    return sum;
#else
    return 1.;
#endif
  }

  smlfloat
  group_select_cats_expect() const {
#if 0
    const unsigned st = paramset_start(PARAMSET_GROUP_SELECT_CATS);
    const unsigned s = paramset_size(PARAMSET_GROUP_SELECT_CATS);

    simple_array<prob_t> pgsc(s);
    typed_cat_pdistro(CAT_PARAM_TYPE::SEL_STRENGTH,CAT_MIX_TYPE::GROUP,pgsc.ptr());

    smlfloat sum(0.);
    for(unsigned i(0);i<(s);++i){ sum += _param[st+i]*pgsc[i]; }
    return sum;
#else
    return 1.;
#endif
  }

  void reset_random_cat_prob(const paramset_t p);

  static
  unsigned
  paramsize_select_model_unscaled(const SELECT_MODEL::index_t sm);

  static
  unsigned
  paramsize_select_model(const SELECT_MODEL::index_t sm);

  unsigned
  paramset_aa_select_matrix_cat_size(const unsigned cat_no) const {
    return paramsize_select_model(select_model(cat_no));
  }

  unsigned
  paramset_aa_select_matrix_cat_start(const unsigned cat_no) const {
    unsigned val(paramset_start(PARAMSET_AA));

    for(unsigned i(0);i<cat_no;++i){
      val += paramset_aa_select_matrix_cat_size(i);
    }

    return val;
  }

  prob_t* bg_pdistro_nscodon_bg_cat(const unsigned bg_cat_no) {
    return _data.bg_pdistro_nscodon_bg_cat[bg_cat_no];
  }

  prob_t* bg_pdistro_nuc_pos_cat(const unsigned cat_no,
                                 const unsigned pos){
    return bg_pdistro_nuc_pos_bg_cat(bg_cat_no_from_cat_no(cat_no),pos);
  }

  prob_t* bg_pdistro_nuc_pos_bg_cat(const unsigned bg_cat_no,
                                    const unsigned pos){
    return _data.bg_pdistro_nuc_pos_bg_cat[bg_cat_no][pos];
  }

  void bg_check() const {
    if(!_data.is_bg_valid) die("rate_gtor_nscodon_base: invalid pdistro bg state");
  }

  // to simplify copy ctor:
  struct auto_copy {
    auto_copy() : is_bg_valid(false) {}

    rate_gtor_nscodon_options_base copt;

    simple_init_matrix<prob_t> bg_pdistro_nscodon_bg_cat;
    simple_init_matrix3d<prob_t> bg_pdistro_nuc_pos_bg_cat;
    simple_init_matrix3d<prob_t> bg_pdistro_nuc_pos_conditioned_on_5p_bg_cat;
    simple_init_matrix3d<prob_t> bg_pdistro_nuc_pos_conditioned_on_3p_bg_cat;

    bool is_bg_valid;
  };

  auto_copy _data;

  // static data:
  //
  struct aa_exchange_param_info {
    aa_exchange_param_info();

    bool is_1nuc_aa_exchange[NSAA::SIZE*NSAA::SIZE];
    unsigned asymm_size;
    unsigned symm_size;
    unsigned asymm_aatoparam_map[NSAA::SIZE*NSAA::SIZE];
    unsigned symm_aatoparam_map[NSAA::SIZE*NSAA::SIZE];
    std::vector<unsigned> asymm_paramtoaa_map;
    std::vector<unsigned> symm_paramtoaa_map;
  };

  struct aa_codon_map_info {
    aa_codon_map_info();

    unsigned nscodon_normal2aa_order_map[NSCODON::SIZE];
    NSCODON::index_t nscodon_aa2normal_order_map[NSCODON::SIZE];
    unsigned aa_codon_size[NSAA::SIZE];
    unsigned aa_codon_start_offset[NSAA::SIZE];
  };

  static const aa_exchange_param_info _aaex;
  static const aa_codon_map_info _aacodon;
};


#endif

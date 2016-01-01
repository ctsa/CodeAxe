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

// $Id: rate_gtor.h 1183 2008-03-27 02:12:11Z ctsa $

/// \file

#ifndef __RATE_GTOR_H
#define __RATE_GTOR_H

#include "observer.h"
#include "param_init_type.h"
#include "param_object.h"
#include "rate_gtor_options.h"
#include "rates_func_options.h"
#include "subs_ml_model_comex.h"
#ifdef DEBUG
#include "subs_ml_ptol.h"
#endif
#include "subs_ml_types.h"
#include "util/general/uncopyable.h"

#include <cstdlib>

#include <iosfwd>

struct cat_manager;
struct condition_func;



/// \brief parameterized rate matrix class
///
struct rate_gtor : public param_object, public observer, private uncopyable {

  typedef param_object base_t;
  typedef rate_gtor self_t;

  rate_gtor(const rate_gtor_options& ropt,
            const rate_gtor_sml_share& smls)
    : base_t(),
      _smls(smls),_is_param_init(false), _ropt(ropt) {
    observe_notifier(_smls.get_bg_gtor_notifier());
  }

  rate_gtor(const self_t& s,
            const rate_gtor_sml_share& smls)
    : base_t(s),
      _smls(smls),_is_param_init(s._is_param_init), _ropt(s._ropt) {
    observe_notifier(_smls.get_bg_gtor_notifier());
  }

  virtual rate_gtor* clone(const rate_gtor_sml_share& smls) const = 0;

  RATE_GTOR_MODEL::index_t model_type() const { return _smls.rate_gtor_model(); }

  unsigned state_size() const { return _smls.state_size(); }

  virtual
  unsigned state_size_conditioned() const { return state_size(); }

  virtual
  unsigned submodel_size() const { return 1; }

  virtual
  unsigned submodel_state_size(unsigned s) const { if(s>=submodel_size()) abort(); return state_size(); }

  virtual
  condition_func* submodel_condition_func_factory(unsigned s) const {
    if(s>=submodel_size()) die("invalid submodel_condition_func_factory");
    return condition_func_factory();
  }

  virtual
  void submodel_pdistro_reduction(const prob_t* in,unsigned s,prob_t* out) const {
    if(s>=submodel_size()) die("invalid pdistro_reduction call");
    std::copy(in,in+state_size(),out);
  }

  virtual
  void submodel_fuse_pdistros(const prob_t* const * in,
                              prob_t* out) const {
    if(submodel_size()>1) die("invalid fuse_pdistros call");
    std::copy(in[0],in[0]+state_size(),out);
  }

  virtual
  void submodel_state_reduction_map(unsigned s,unsigned* out) const {
    if(s>=submodel_size()) abort();
    const unsigned ss(state_size());
    for(unsigned i(0);i<ss;++i) out[i] = i;
    // remap ambiguous state too:
    out[ss] = ss;
  }

  virtual
  void submodel_adjust_p(prob_t*,const unsigned) const {}

  SITE_MODEL::index_t site_model() const {
    const SITE_MODEL::index_t sm(RATE_GTOR_MODEL::convert_to_site_model(model_type()));
    if(sm==SITE_MODEL::NONE) die("invalid site-model");
    return sm;
  }

  virtual void store_state(std::ostream&) const {}

  void load_state(std::istream& is);

  void
  dump(const char* label,
       std::ostream& os) const;

  virtual
  void
  report(std::ostream& os) const = 0;

  virtual
  condition_func* condition_func_factory() const;

  /// \brief get the instantaneous rates
  virtual void rates(smlfloat rates[],
                     const rates_func_options& opt) const = 0;

  void reset_param(const PARAM_INIT_TYPE::index_t pinit){
    reset_param_internal(pinit);
    _is_param_init=true;
    fix_input_param();
  }

  // per site lnp penalty
  virtual
  smlfloat param_penalty() const { return 0.;}

  const prob_t* bg_pdistro_cat(const unsigned cat_no) const {
    return bg_pdistro_bg_cat(bg_cat_no_from_cat_no(cat_no));
  }

  // fix_input_param sets input parameters to legal values (normalizing
  // prob distros -- it is intended to repair input parameters and will
  // change the meaning of the model.
  //
  void fix_input_param();

  /////////////////////////////////////////////////////////////////
protected:

  virtual
  void reset_param_internal(const PARAM_INIT_TYPE::index_t pinit);

  virtual void bg_pdistro_update_internal() {}

  virtual void load_state_internal(std::istream&) {}

  virtual void fix_input_param_internal() {}

  virtual void param_update_from_new_bg_pdistro_internal() {}

  unsigned bg_cat_size() const {
    return _smls.bg_cat_size();
  }

  unsigned bg_cat_no_from_cat_no(const unsigned cat_no) const {
    return _smls.bg_cat_no_from_cat_no(cat_no);
  }

  const prob_t* bg_pdistro_bg_cat(const unsigned bg_cat_no) const {
    return _smls.bg_state_pdistro_bg_cat(bg_cat_no);
  }

  unsigned
  param_end() const { return 0; }

  const cat_manager& cm() const { return _smls.get_cat_manager(); }

private:
  //////////////////////////////////////////////////////////
  template <typename FloatType>
  void
  set_uniform_param_state(FloatType x) {
    set_param_state(x);
    fix_input_param();
  }

  void check_param() const;

  void param_update_from_new_bg_pdistro() {
    if(_is_param_init){
      param_update_from_new_bg_pdistro_internal();
      check_param();
    }
  }

  void bg_pdistro_update() {
    bg_pdistro_update_internal();
    param_update_from_new_bg_pdistro();
  }

  virtual void recieve_notifier_event(const notifier* n,
                                      const EVENT_TYPE::index_t e) {
    if(n == &_smls.get_bg_gtor_notifier() && e == EVENT_TYPE::PDISTRO_CHANGE){
      bg_pdistro_update();
    }
  }


  ///////////////////////////
  const rate_gtor_sml_share _smls; ///< store external reference to shared subs_ml_model data
  bool _is_param_init;
  rate_gtor_options _ropt;
};

#endif

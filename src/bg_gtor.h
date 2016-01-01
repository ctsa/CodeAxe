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

// $Id: bg_gtor.h 1153 2008-03-18 00:24:07Z ctsa $

/// \file

#ifndef __BG_GTOR_H
#define __BG_GTOR_H

#include "bg_gtor_options.h"
#include "observer.h"
#include "param_composite.h"
#include "param_init_type.h"
#include "prob_gtor_param.h"
#include "prob_gtor_obs.h"
#include "rate_gtor_options.h"
#include "subs_ml_model_comex.h"
#include "subs_ml_types.h"
#include "util/general/uncopyable.h"

#include <cassert>

#include <iosfwd>
#include <vector>

struct cat_manager;



/// \brief subs_ml_model rate bg distribution handler
///
/// distro_cats = the actual number of different distros, distro cats is the only
/// number external clients should worry about, or they can ask for conversions
/// from the global cat no.
///
/// bg_cats = the number of distro cats if there are bg params, otherwise
/// this number is irrelevent. client code shouldn't have access to this number
///
struct bg_gtor : public param_composite, public observer, public notifier, private uncopyable {

  typedef param_composite base_t;
  typedef bg_gtor self_t;

  typedef prob_gtor ptor_t;
  typedef prob_gtor_param ptorp_t;
  typedef prob_gtor_obs ptoro_t;

  bg_gtor(std::istream& is,
          const bg_gtor_sml_share& smls)
    : _smls(smls) {
    load_state(is);
    register_param();
    observe_notifier(_smls.get_obs_info_notifier());
  }

  bg_gtor(const bg_gtor_options bgo,
          const bg_gtor_sml_share& smls)
    : _smls(smls) {
    _data.bgo=bgo;
    init();
    register_param();
    observe_notifier(_smls.get_obs_info_notifier());
  }

  bg_gtor(const self_t& s,
          const bg_gtor_sml_share& smls);

  ~bg_gtor() {
    const unsigned ds(distro_cat_size());
    for(unsigned i(1);i<=ds;++i) delete _ptor[ds-i];
  }

private:
  enum ptor_type {
    PT_PARAM,
    PT_PARAM_LINKED,
    PT_OBS
  };

public:

  unsigned
  state_size() const { return _smls.state_size(); }

  /// set all parameterized_distro_cats as close as possible to x:
  ///
  /// intentionally adjusts all params, including no-train
  template <typename ForwardIterator>
  void set_pdistro(ForwardIterator x) {
    if(ptype()==PT_OBS) return;
    const unsigned ds(distro_cat_size());
    for(unsigned i(0);i<ds;++i){ set_distro_cat_state_pdistro(i,x); }
    notify_observers(EVENT_TYPE::PDISTRO_CHANGE);
  }

private:
  /// set parameterized_distro_cat as close as possible to x
  ///
  /// intentionally adjusts all params, including no-train
  template <typename InputIterator>
  void set_distro_cat_state_pdistro(const unsigned distro_cat_no,
                                    InputIterator x) {
    if(ptype()==PT_OBS) return;

    // for param roots, distro_cat_no == root_cat_no:
    distro_cat_ptorp(distro_cat_no).set_pdistro(x);
  }
public:

  unsigned distro_cat_size() const { return _data.distro_cat_size; }

  const prob_t*
  distro_cat_state_pdistro(const unsigned distro_cat_no) const {
    return static_cast<const ptor_t*>(_ptor[distro_cat_no])->pdistro();
  }

  unsigned
  get_distro_cat_no_from_cat_no(const unsigned cat_no) const {
    return _data.distro_cat_map[cat_no];
  }

  const prob_t*
  cat_state_pdistro(const unsigned cat_no) const {
    return distro_cat_state_pdistro(get_distro_cat_no_from_cat_no(cat_no));
  }

  void reset_param(const PARAM_INIT_TYPE::index_t pinit);

  void store_state(std::ostream& os) const;

private:
  void load_state(std::istream& is);

  ptorp_t&
  distro_cat_ptorp(const unsigned distro_cat_no) {
    return static_cast<ptorp_t&>(*_ptor[distro_cat_no]);
  }

  const ptorp_t&
  distro_cat_ptorp(const unsigned distro_cat_no) const {
    return static_cast<const ptorp_t&>(*_ptor[distro_cat_no]);
  }

  ptoro_t&
  distro_cat_ptoro(const unsigned distro_cat_no) {
    assert(distro_cat_no<distro_cat_size());
    return static_cast<ptoro_t&>(*_ptor[distro_cat_no]);
  }

  const ptoro_t&
  distro_cat_ptoro(const unsigned distro_cat_no) const {
    assert(distro_cat_no<distro_cat_size());
    return static_cast<const ptoro_t&>(*_ptor[distro_cat_no]);
  }

  void
  init() {
    setup_distro_cats();
    const unsigned ds(distro_cat_size());
    _ptor.resize(ds);
    if(ptype()==PT_OBS){
      for(unsigned i(0);i<ds;++i) _ptor[i]=prob_gtor_obs_factory(i);
    } else {
      for(unsigned i(0);i<ds;++i) _ptor[i]=prob_gtor_param_factory();
    }
  }

  void
  register_param() {
    clear_register();
    if(ptype()==PT_OBS) return;
    const unsigned ds(distro_cat_size());
    for(unsigned i(0);i<ds;++i){
      base_t::register_param(distro_cat_ptorp(i));
    }
  }

  const cat_manager& cm() const {
    return _smls.get_cat_manager();
  }

  void setup_distro_cats();

  ptor_type
  ptype() const {
    using namespace SUBS_RATE_BG_MODEL;
    switch(_data.bgo.bg){
    case OBS_AVG:        return PT_OBS;
    default:             return PT_PARAM;
    }
  }

  ptorp_t* prob_gtor_param_factory() const;

  ptoro_t* prob_gtor_obs_factory(const unsigned distro_cat_no) const;

  virtual
  void recieve_notifier_event(const notifier* n,
                              const EVENT_TYPE::index_t e){
    if(n == &_smls.get_obs_info_notifier() && e == EVENT_TYPE::PDISTRO_CHANGE){
      if(ptype()==PT_OBS){
        notify_observers(EVENT_TYPE::PDISTRO_CHANGE);
      }
    }
  }

  virtual
  void
  post_write_param_state() {
    if(ptype()==PT_PARAM){
      notify_observers(EVENT_TYPE::PDISTRO_CHANGE);
    }
  }

  struct auto_copy {
    bg_gtor_options bgo;
    std::vector<unsigned> distro_cat_map;
    unsigned distro_cat_size;
  };

  const bg_gtor_sml_share _smls;
  std::vector<ptor_t*> _ptor;
  auto_copy _data;
};

#endif

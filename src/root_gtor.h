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

// $Id: root_gtor.h 1123 2008-01-29 00:35:28Z ctsa $

/// \file

#ifndef __ROOT_GTOR_H
#define __ROOT_GTOR_H

#include "param_composite.h"
#include "param_init_type.h"
#include "prob_gtor_param.h"
#include "prob_gtor_obs.h"
#include "rate_gtor_options.h"
#include "root_gtor_options.h"
#include "subs_ml_model_comex.h"
#include "subs_ml_types.h"
#include "util/general/uncopyable.h"

#include <cassert>

#include <iosfwd>
#include <vector>

struct cat_manager;


/// \brief subs_ml_model root distribution handler
///
/// note that the number of "distro categories" will be >= to the
/// number of root categories. If each root cat is parameterized,
/// root cats == distro cats. When a root cat is set to be derived
/// from other parameters in the model, each root cat could have a
/// number of distro cats of up to the the number of different global
/// categories associated with that root cat. from the users
/// perspective all distros must be requested by the global cat
/// number.
///
/// distro_cats = the actual number of different distros, distro cats is the only
/// number external clients should worry about, or they can ask for conversions
/// from the global cat no.
///
/// root_cats = the number of distro cats if there are root params, otherwise
/// this number is irrelevent. client code shouldn't have access to this number
///
struct root_gtor : public param_composite, private uncopyable {

  typedef param_composite base_t;
  typedef root_gtor self_t;

  typedef prob_gtor ptor_t;
  typedef prob_gtor_param ptorp_t;
  typedef prob_gtor_obs ptoro_t;

  root_gtor(std::istream& is,
            const root_gtor_sml_share& smls)
    : _smls(smls) {
    load_state(is);
    register_param();
  }

  root_gtor(const ROOT_GTOR_MODEL::index_t rogm,
            const root_gtor_sml_share& smls)
    : _smls(smls) {
    _data.rogm=rogm;
    init();
    register_param();
  }

  root_gtor(const self_t& s,
            const root_gtor_sml_share& smls);

  ~root_gtor() {
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

  /// set all parameterized_distro_cats/root_cats as close as possible to x:
  ///
  /// intentionally adjusts all params, including no-train
  template <typename ForwardIterator>
  void set_pdistro(ForwardIterator x) {
    if(ptype()==PT_OBS) return;
    const unsigned ds(distro_cat_size());
    for(unsigned i(0);i<ds;++i){ set_distro_cat_state_pdistro(i,x); }
  }

  /// set parameterized_distro_cat/root cat as close as possible to x
  ///
  /// intentionally adjusts all params, including no-train
  template <typename InputIterator>
  void set_distro_cat_state_pdistro(const unsigned distro_cat_no,
                                    InputIterator x) {
    if(ptype()==PT_OBS) return;

    // for param roots, distro_cat_no == root_cat_no:
    distro_cat_ptorp(distro_cat_no).set_pdistro(x);
  }

  unsigned distro_cat_size() const { return _data.distro_cat_size; }

  const prob_t*
  distro_cat_state_pdistro(const unsigned distro_cat_no) const {
    return static_cast<const ptor_t*>(_ptor[distro_cat_no])->pdistro();
  }

  void
  report(const char* label,
         std::ostream& os) const;

  unsigned
  get_distro_cat_no_from_cat_no(const unsigned cat_no) const {
    return _data.distro_cat_map[cat_no];
  }

  const prob_t*
  cat_state_pdistro(const unsigned cat_no) const {
    return distro_cat_state_pdistro(get_distro_cat_no_from_cat_no(cat_no));
  }

  /// these functions enable the rootcycle minimizer
  ///
  param_object_base&
  cat_param_object(const unsigned cat_no) {
    if(ptype()==PT_OBS) return *this;
    else                return distro_cat_ptorp(get_distro_cat_no_from_cat_no(cat_no));
  }

  const param_object_base&
  cat_param_object(const unsigned cat_no) const {
    if(ptype()==PT_OBS) return *this;
    else                return distro_cat_ptorp(get_distro_cat_no_from_cat_no(cat_no));
  }

  void reset_param(const PARAM_INIT_TYPE::index_t pinit);

  void store_state(std::ostream& os) const;

  // poorly named -- notifies that the enclosing model instance
  // (subs_ml_model) has just completed a parameter update.
  //
  /// \todo distinguish between parent parameter update and obs_info update
  ///
  void parent_model_update() {
    if(ptype()!=PT_OBS) return;
    const unsigned ds(distro_cat_size());
    for(unsigned i(0);i<ds;++i) distro_cat_ptoro(i).pdistro_update();
  }

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
      for(unsigned i(0);i<ds;++i) _ptor[i]=prob_gtor_param_factory(i);
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
    using namespace ROOT_GTOR_MODEL;
    switch(_data.rogm){
    case OBS_AVG:
    case OBS_TIME_AVG:
    case LSPROB:        return PT_OBS;
    case FULL_GC_SHIFT: return PT_PARAM_LINKED;
    default:            return PT_PARAM;
    }
  }

  ptorp_t* prob_gtor_param_factory(const unsigned distro_cat_no) const;

  ptoro_t* prob_gtor_obs_factory(const unsigned distro_cat_no) const;

  struct auto_copy {
    ROOT_GTOR_MODEL::index_t rogm;
    std::vector<unsigned> distro_cat_map;
    unsigned distro_cat_size;
  };

  const root_gtor_sml_share _smls;
  std::vector<ptor_t*> _ptor;
  auto_copy _data;
};

#endif

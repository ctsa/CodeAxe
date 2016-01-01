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

// $Id: obs_info.h 1121 2008-01-28 22:32:57Z ctsa $

/// \file

#ifndef __OBS_INFO_H
#define __OBS_INFO_H

#include "observer.h"
#include "param_init_type.h"
#include "simple_util.h"
#include "subs_ml_types.h"
#include "util/general/die.h"

#include <algorithm>
#include <iosfwd>

struct obs_info : public notifier {

  enum source_type {
    NONE,
    COUNTS,
    DISTRO
  };

  obs_info(const unsigned n_obs_cats,
           const unsigned n_seqs,
           const unsigned n_states)
    : _source(NONE), _data(n_obs_cats,n_seqs,n_states) {}

  obs_info(const unsigned n_obs_cats,
           const unsigned n_seqs,
           const unsigned n_states,
           std::istream& is)
    : _source(NONE), _data(n_obs_cats,n_seqs,n_states) { load_state(is); }

  void store_state(std::ostream& os) const;

  void
  set_obs_cat_seq_state_distro(prob_t const * const * const * c);

  void
  set_obs_cat_seq_state_counts(smlfloat const * const * const * c);

  const prob_t * const * seq_state_distro() const {
    ucheck();
    return _data.seq_state_distro.ptr();
  }

  const prob_t * state_distro() const {
    ucheck();
    return _data.state_distro.ptr();
  }

  const prob_t * const * const * seq_state_distro_obs_cat() const {
    ucheck();
    return _data.seq_state_distro_obs_cat.ptr();
  }

  const prob_t * const * state_distro_obs_cat() const {
    ucheck();
    return _data.state_distro_obs_cat.ptr();
  }

  void reset_param(const PARAM_INIT_TYPE::index_t pinit);

private:
  template <typename T>
  void
  set_obs_cat_seq_state_int(T const * const * const * c);

  void update_distros();

  void load_state(std::istream& is);

  void
  ucheck() const {
    if(_source==NONE) die("Invalid obs_info state");
  }

  /////////////////////////////////////////
  struct obs_data {

    obs_data(const unsigned a,
             const unsigned b,
             const unsigned c)
      : cat_seq_state_source(a,b,c),
        seq_state_distro(b,c),
        state_distro(c),
        seq_state_distro_obs_cat(a,b,c),
        state_distro_obs_cat(a,c){}

    simple_matrix3d<smlfloat> cat_seq_state_source; // this can be a distro or (fractional) counts
    simple_matrix<prob_t> seq_state_distro;
    simple_array<prob_t> state_distro;
    simple_matrix3d<prob_t> seq_state_distro_obs_cat;
    simple_matrix<prob_t> state_distro_obs_cat;
  };

  source_type _source;
  obs_data _data;
};


#endif

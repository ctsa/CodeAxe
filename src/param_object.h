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

// $Id: param_object.h 1113 2008-01-27 00:23:59Z ctsa $

/// \file

#ifndef __PARAM_OBJECT_H
#define __PARAM_OBJECT_H

#include "param_object_base.h"
#include "util/general/die.h"

#include <vector>



/// \todo -- clarify the non-negative meaning -- as is, things get mirrored if
/// values are coming from a minimizer, but an error occurs if values are coming from a human, right?
///
///
struct param_object : public param_object_base {

  typedef param_object_base base_t;
  typedef param_object self_t;

  explicit
  param_object(const unsigned init_param_size=0)
    : base_t(),_param(init_param_size,0),
      _is_train_param(init_param_size,true),
      _is_free_block_start(init_param_size,false),
      _is_free_block_stop(init_param_size,false),
      _is_nonnegative(init_param_size,true),
      _param_variance() {}

  param_object(const self_t& s)
    : base_t(s), _param(s._param),
      _is_train_param(s._is_train_param),
      _is_free_block_start(s._is_free_block_start),
      _is_free_block_stop(s._is_free_block_stop),
      _is_nonnegative(s._is_nonnegative),
      _param_variance() {}

  virtual
  unsigned param_size(const PARAM_VIEW::index_t pv = PARAM_VIEW::ALL) const;

  virtual
  smlfloat*
  param_state(smlfloat* x,
              const PARAM_VIEW::index_t pv = PARAM_VIEW::ALL) const;

  virtual
  const smlfloat*
  set_param_state(const smlfloat* x,
                  const PARAM_VIEW::index_t pv = PARAM_VIEW::ALL);

  virtual
  void
  set_param_state(const smlfloat x);

  virtual
  bool*
  is_train_param_state(bool* x) const {
    const unsigned ps(param_size());
    for(unsigned i(0);i<ps;++i) {
      *x = _is_train_param[i];
      ++x;
    }
    return x;
  }

  virtual
  const bool*
  set_is_train_param_state(const bool* x) {
    const unsigned ps(param_size());
    for(unsigned i(0);i<ps;++i) {
      _is_train_param[i] = *x;
      ++x;
    }
    return x;
  }

  virtual
  void
  set_is_train_param_state(const bool b) {
    const unsigned ps(param_size());
    for(unsigned i(0);i<ps;++i) _is_train_param[i] = b;
  }

  virtual
  const smlfloat*
  attach_ci(const smlfloat* ci){
    const unsigned ps(param_size());

    _param_variance.clear();
    _param_variance.resize(ps,0.);

    for(unsigned i(0);i<ps;++i){
      if(_is_train_param[i]){
        _param_variance[i] = *ci;
        ci++;
      }
    }
    return ci;
  }

  // public unit functions:
  //
  /// \todo add more of these!
  virtual
  smlfloat
  param_variance(unsigned i) const {
    if(_param_variance.size()>i){
      return _param_variance[i];
    } else {
      return smlfloat(0.);
    }
  }


protected:
  virtual
  unsigned
  pseudo_param_size() const { return 0; }

  virtual
  smlfloat
  trainable_param_transform(const unsigned i) const {
    return _param[i];
  }

  virtual
  void
  trainable_param_untransform(const unsigned i,
                              const smlfloat p) {
    _param[i] = p;
  }

  virtual
  void
  set_pseudo_param(const unsigned param_no,
                   const smlfloat /*val*/){
    if(param_no>=pseudo_param_size()){
      die("Bad pseudoparam number");
    }
  }

  void
  resize(const unsigned p){
    _param.resize(p,0);
    _is_train_param.resize(p,true);
    _is_free_block_start.resize(p,false);
    _is_free_block_stop.resize(p,false);
    _is_nonnegative.resize(p,true);
  }

protected:
  virtual
  void
  post_write_param_state() {}

protected:
  std::vector<smlfloat> _param;
  std::vector<bool> _is_train_param;
  std::vector<bool> _is_free_block_start;
  std::vector<bool> _is_free_block_stop;
  std::vector<bool> _is_nonnegative;
private:
  std::vector<smlfloat> _param_variance;
};

#endif

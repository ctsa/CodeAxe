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

// $Id: param_composite.h 1153 2008-03-18 00:24:07Z ctsa $

/// \file

#ifndef __PARAM_COMPOSITE_H
#define __PARAM_COMPOSITE_H

#include "param_object_base.h"
#include "util/general/die.h"

#include <vector>



///
struct param_composite : public param_object_base {

  typedef param_object_base base_t;
  typedef param_composite self_t;

  param_composite() : base_t() {}

  param_composite(const self_t& s) : base_t(s) {}

  virtual
  unsigned
  param_size(const PARAM_VIEW::index_t pv = PARAM_VIEW::ALL) const;

  virtual
  smlfloat*
  param_state(smlfloat* x,
              const PARAM_VIEW::index_t pv = PARAM_VIEW::ALL) const {

    const unsigned s(_obj.size());
    for(unsigned i(0);i<s;++i){
      x = _obj[i]->param_state(x,pv);
    }
    return x;
  }

  virtual
  const smlfloat*
  set_param_state(const smlfloat* x,
                  const PARAM_VIEW::index_t pv = PARAM_VIEW::ALL) {

    const unsigned s(_obj.size());
    for(unsigned i(0);i<s;++i){
      x = _obj[i]->set_param_state(x,pv);
    }
    post_write_param_state();
    return x;
  }

  virtual
  void
  set_param_state(const smlfloat x) {

    const unsigned s(_obj.size());
    for(unsigned i(0);i<s;++i){
      _obj[i]->set_param_state(x);
    }
    post_write_param_state();
  }

  virtual
  bool*
  is_train_param_state(bool* x) const {

    const unsigned s(_obj.size());
    for(unsigned i(0);i<s;++i){
      x = _obj[i]->is_train_param_state(x);
    }
    return x;
  }

  virtual
  const bool*
  set_is_train_param_state(const bool* x) {

    const unsigned s(_obj.size());
    for(unsigned i(0);i<s;++i){
      x = _obj[i]->set_is_train_param_state(x);
    }
    return x;
  }

  virtual
  void
  set_is_train_param_state(const bool b);

  virtual
  const smlfloat*
  attach_ci(const smlfloat* x) {

    const unsigned s(_obj.size());
    for(unsigned i(0);i<s;++i){
      x = _obj[i]->attach_ci(x);
    }
    return x;
  }

  virtual
  smlfloat
  param_variance(unsigned n) const {
    const unsigned s(_obj.size());
    for(unsigned i(0);i<s;++i){
      const unsigned ps(_obj[i]->param_size());
      if(n>=ps) n -= ps;
      else      return _obj[i]->param_variance(n);
    }
    die("param_composite: Invalid param_variance call");
  }

protected:
  void
  clear_register(){
    _obj.clear();
  }

  void
  register_param(param_object_base& o){
    _obj.push_back(&o);
  }

private:
  virtual
  void
  post_write_param_state()  {}

  std::vector<param_object_base*> _obj;
};


#endif

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

// $Id: prob_gtor_param.h 1120 2008-01-28 21:24:46Z ctsa $

/// \file

#ifndef __PROB_GTOR_PARAM_H
#define __PROB_GTOR_PARAM_H

#include "observer.h"
#include "prob_gtor.h"
#include "param_object.h"
#include "util/math/prob_util.h"

#include <iosfwd>
#include <set>


/// \brief
///
struct prob_gtor_param :  public prob_gtor, public param_object, public notifier {

  typedef param_object obase_t;
  typedef prob_gtor pbase_t;
  typedef notifier nbase_t;
  typedef prob_gtor_param self_t;

  prob_gtor_param(const unsigned init_state_size,
                  const unsigned init_param_size)
    : pbase_t(init_state_size), obase_t(init_param_size), nbase_t() {
    if(_is_free_block_start.size()) _is_free_block_start.front() = true;
    if(_is_free_block_stop.size()) _is_free_block_stop.back() = true;
  }

  prob_gtor_param(const self_t& s) : pbase_t(s), obase_t(s), nbase_t(s) {}

  virtual self_t* clone() const = 0;

  template <typename InputIterator>
  InputIterator
  set_pdistro(InputIterator x) {
    const unsigned ss(state_size());
    for(unsigned i(0);i<ss;++i) { pdistro()[i]=*x; ++x;}

    // recycle distro into distro set from parameters:
    get_param_from_pdistro();
    update_from_new_param();;
    return x;
  }

  void
  store_state(const char* label,
              std::ostream& os) const;

private:
  void
  update_from_new_param(){
    param_norm();
    get_pdistro_from_param();
    notify_observers(EVENT_TYPE::PDISTRO_CHANGE);
  }

  virtual
  void get_param_from_pdistro() = 0;

  virtual
  void param_norm() {
    pdistro_norm_free_param(_param.begin(),
                            _param.end(),
                            _is_train_param.begin());
  }

  virtual
  void get_pdistro_from_param() = 0;

  virtual
  void
  post_write_param_state() {
    update_from_new_param();
  }
};


#endif

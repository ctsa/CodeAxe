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

// $Id: param_object_base.h 1114 2008-01-27 01:13:43Z ctsa $

/// \file

#ifndef __PARAM_OBJECT_BASE_H
#define __PARAM_OBJECT_BASE_H

#include "param_view.h"
#include "subs_ml_types.h"


struct param_object_base {

  typedef param_object_base self_t;

  virtual ~param_object_base() {}

private:
  self_t& operator=(const self_t&);
public:

  virtual
  unsigned
  param_size(const PARAM_VIEW::index_t pv = PARAM_VIEW::ALL) const = 0;

  virtual
  smlfloat*
  param_state(smlfloat* x,
              const PARAM_VIEW::index_t pv = PARAM_VIEW::ALL) const = 0;

  virtual
  const smlfloat*
  set_param_state(const smlfloat* x,
                  const PARAM_VIEW::index_t pv = PARAM_VIEW::ALL) = 0;

  virtual
  void
  set_param_state(const smlfloat x) = 0;

  virtual
  bool*
  is_train_param_state(bool* x) const = 0;

  virtual
  const bool*
  set_is_train_param_state(const bool* x) = 0;

  virtual
  void
  set_is_train_param_state(const bool b) = 0;

  virtual
  const smlfloat*
  attach_ci(const smlfloat* x) = 0;

  virtual
  smlfloat
  param_variance(unsigned n) const = 0;
};


#endif

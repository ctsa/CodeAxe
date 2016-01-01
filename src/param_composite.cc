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

// $Id: param_composite.cc 1111 2008-01-26 01:11:52Z ctsa $

/// \file

#include "param_composite.h"


unsigned
param_composite::
param_size(const PARAM_VIEW::index_t pv) const {

  unsigned ps(0);
  const unsigned os(_obj.size());
  for(unsigned i(0);i<os;++i) ps += _obj[i]->param_size(pv);
  return ps;
}



void
param_composite::
set_is_train_param_state(const bool b) {

  const unsigned s(_obj.size());
  for(unsigned i(0);i<s;++i){
    _obj[i]->set_is_train_param_state(b);
  }
}

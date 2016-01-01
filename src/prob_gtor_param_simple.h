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

// $Id: prob_gtor_param_simple.h 1078 2008-01-07 02:37:41Z ctsa $

/// \file

#ifndef __PROB_GTOR_PARAM_SIMPLE_H
#define __PROB_GTOR_PARAM_SIMPLE_H


#include "prob_gtor_param.h"


/// \brief simple prob generator: 1 parameter per distro element
///
struct prob_gtor_param_simple : public prob_gtor_param {

  typedef prob_gtor_param_simple self_t;
  typedef prob_gtor_param base_t;

  explicit
  prob_gtor_param_simple(const unsigned init_state_size)
    : base_t(init_state_size,init_state_size) {}

  virtual self_t* clone() const { return new self_t(*this); }

private:
  virtual
  void get_param_from_pdistro() {
    const unsigned ps(param_size());
    for(unsigned i(0);i<ps;++i) _param[i]= pdistro()[i];
  }

  virtual
  void get_pdistro_from_param() {
    const unsigned ps(param_size());
    for(unsigned i(0);i<ps;++i) pdistro()[i] = _param[i];
  }
};


#endif

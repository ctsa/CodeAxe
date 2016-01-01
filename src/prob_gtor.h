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

// $Id: prob_gtor.h 1113 2008-01-27 00:23:59Z ctsa $

/// \file

#ifndef __PROB_GTOR_H
#define __PROB_GTOR_H

#include "simple_util.h"
#include "subs_ml_types.h"
#include "util/general/die.h"


/// \brief probability generator base class
///
struct prob_gtor {

  typedef prob_gtor self_t;

  explicit
  prob_gtor(const unsigned init_state_size) : _distro(init_state_size) {}

  virtual ~prob_gtor() {}

private:
  self_t& operator=(const self_t&);

public:

  virtual
  const prob_t* pdistro() const { return _distro.ptr(); }

  unsigned state_size() const { return _distro.size(); }

protected:

  prob_t* pdistro() { return _distro.ptr();}
  prob_t* pdistro_mutable() const { return _distro.ptr();}

private:
  mutable simple_array<prob_t> _distro;
};


#endif

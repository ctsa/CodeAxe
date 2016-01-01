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

// $Id: estep_dat_t.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __ESTEP_DAT_T_H
#define __ESTEP_DAT_T_H

#include "subs_ml_types.h"

/// \brief data calculated and stored in the em e-step
///   trans = matrix of expected substitutions per branch
///   root = (posterior prob of branch node state?)
///
struct estep_dat_t {

  estep_dat_t(const unsigned init_n_states,
              const unsigned init_n_branches) :
    n_states(init_n_states),
    n_branches(init_n_branches),
    root(new smlfloat[n_states]),
    trans(new smlfloat*[n_branches]){
    const unsigned n2(n_states*n_states);

    trans[0] = new smlfloat[n_branches*n2];
    for(unsigned i(0);i<n_branches;++i){ trans[i] = trans[0]+i*n2; }
  }

  ~estep_dat_t() {
    delete [] root;
    delete [] trans[0];
    delete [] trans;
  }

  unsigned n_states;
  unsigned n_branches;
  smlfloat* root;
  smlfloat** trans;
};


#endif

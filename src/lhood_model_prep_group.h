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

// $Id: lhood_model_prep_submodel.h 1050 2007-12-05 05:34:52Z ctsa $

/// \file

#ifndef __LHOOD_MODEL_PREP_GROUP_H
#define __LHOOD_MODEL_PREP_GROUP_H

#include "subs_ml_types.h"
#include "util/general/uncopyable.h"

#include <vector>


namespace PVL {
extern "C" {
#include "pr.h"
}
}



struct group_data : private uncopyable {
  explicit
  group_data(const unsigned n_groups)
    : group_prob(n_groups), group_prob_tmp(n_groups) {}

  std::vector<PVL::pr_t> group_prob;
  std::vector<smlfloat> group_prob_tmp;
};


#endif

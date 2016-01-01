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

// $Id: condition_func.h 1054 2007-12-07 19:37:22Z ctsa $

/// \file

#ifndef __CONDITION_FUNC_H
#define __CONDITION_FUNC_H

#include "subs_ml_types.h"
#include "util/general/uncopyable.h"

#include <vector>

struct bi_tree;
struct site_data_fastlup_core;
struct tree_probs;


struct condition_func : private uncopyable {

  virtual ~condition_func() {}

  /// init condition function for data evaluated by the model
  ///
  virtual
  void
  data_init(const bi_tree&,
            const site_data_fastlup_core&) {}

  /// condition site probs based on the input branch props.
  ///
  void
  condition_site_prob(const tree_probs& tprob,
                      const std::vector<bool>& site_prob_mask,
                      prob_t* site_prob) const {
    tprob_init(tprob);
    condition_site_prob_prep(site_prob_mask,site_prob);
  }

  /// initialize conditioning levels as a function of all model parameters
  ///
  virtual
  void
  tprob_init(const tree_probs&) const {}

  /// update conditioning levels for new root distro.
  ///
  virtual
  void
  root_update(const prob_t*) const {}

  /// condition site probs based on last update of conditioning levels
  ///
  virtual
  void
  condition_site_prob_prep(const std::vector<bool>&,
                           prob_t*) const {}

  virtual
  smlfloat
  transform_prob(const smlfloat p) const { return p; }

};

#endif

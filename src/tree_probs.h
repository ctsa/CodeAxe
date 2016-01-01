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

// $Id: tree_probs.h 1161 2008-03-20 19:28:26Z ctsa $

/// \file

#ifndef __TREE_PROBS_H
#define __TREE_PROBS_H

#include "simple_util.h"
#include "subs_ml_types.h"
#include "util/general/workspace.h"


struct subs_ml_model;


/// \brief convenience struct to hold branch and root probabilities
/// for a tree
///
struct tree_probs {

  typedef tree_probs self_t;

  tree_probs()
    : n_states(0), n_branches(0), _n_sites(0), _is_site_prob_mode(false) {}

  /// initialize probabilities. **matrices are stored transposed**
  /// so that the 'from' state is the fast index: val = m[from+SIZE*to]
  ///
  void
  init_probs(const subs_ml_model& mdl,
             const unsigned cat = 0,
             const bool is_use_submodel = false,
             const unsigned submodel_no = 0) {
    workspace<char> ws[3];
    init_probs(mdl,cat,is_use_submodel,submodel_no,ws);
  }

  /// malloc free version
  void
  init_probs(const subs_ml_model& mdl,
             const unsigned cat,
             const bool is_use_submodel,
             const unsigned submodel_no,
             workspace<char> ws[3]);

  // update only the root distro from the value calculated in init_probs()
  //
  void
  update_root(const prob_t* root);

  void self_check() const;

  prob_t* root_prob() { return _root_prob.ptr(); }
  prob_t const * root_prob() const { return _root_prob.ptr(); }

  prob_t* branch_prob(const unsigned b) { return _branch_prob[b]; }
  prob_t const * branch_prob(const unsigned b) const { return _branch_prob[b]; }

  bool& is_branch_tprob_identity(const unsigned b){
    return _is_branch_tprob_identity[b];
  }

  bool is_branch_tprob_identity(const unsigned b) const {
    return _is_branch_tprob_identity[b];
  }

private:
  void alloc(const unsigned _n_states,
             const unsigned _n_branches);

public:
  unsigned n_states;
  unsigned n_branches;
private:
  unsigned _n_sites;
  bool _is_site_prob_mode;
  simple_init_array<prob_t> _root_prob;
  simple_init_matrix<prob_t> _branch_prob;
  simple_init_array<bool> _is_branch_tprob_identity;
};


#endif

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

// $Id: condition_func_overlap_base.h 1054 2007-12-07 19:37:22Z ctsa $

/// \file

#ifndef __CONDITION_FUNC_OVERLAP_BASE_H
#define __CONDITION_FUNC_OVERLAP_BASE_H


#include "condition_func.h"
#include "tree_probs.h"
#include "simple_util.h"
#include "site_data_fastlup_core.h"

#include <vector>

/// \brief condition function for overlapping sites:
///
struct condition_func_overlap_base : public condition_func {

  condition_func_overlap_base() : tree_ptr(0) {}

  virtual
  void
  data_init(const bi_tree& tree,
            const site_data_fastlup_core& sdf);

  virtual
  void
  tprob_init(const tree_probs& tprob) const;

  virtual
  void
  root_update(const prob_t* root) const;

  virtual
  void
  condition_site_prob_prep(const std::vector<bool>& site_prob_mask,
                           prob_t* site_prob) const;

protected:
  virtual
  unsigned state_size() const {
    return conditioning_state_size()*conditioned_state_size();
  }

private:
  virtual bool is_bidirectional() const { return false; }

  virtual
  void get_leaf_leg(prob_t* legp,
                    const prob_t* tprob,
                    const unsigned conditioning_state,
                    const unsigned condition_no) const = 0;

  virtual unsigned conditioning_state_size() const = 0;
  virtual unsigned conditioned_state_size() const = 0;

  virtual
  unsigned state_2_conditioning_state(unsigned i,
                                      const unsigned condition_no) const = 0;

  struct cf_data {
    mutable simple_init_array<prob_t> condition_prob;
    mutable tree_probs tprob_cf;
    site_data_fastlup_core sdf_cf;
    std::vector<unsigned> data_reduced_state_index_map;
  };

  void
  data_init_cfd(cf_data& cfd,
                const unsigned condition_no,
                const bi_tree& tree,
                const site_data_fastlup_core& sdf);

  void
  update_condition_prob(const cf_data& cfdx) const;

  void
  tprob_init_cfd(const cf_data& cfdx,
                 const unsigned condition_no,
                 const tree_probs& tprob) const;

  void
  root_update_cfd(const cf_data& cfdx,
                  const prob_t* root) const;

  const bi_tree* tree_ptr;
  cf_data cfd;
  cf_data cfd2;
};

#endif

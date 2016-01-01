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

// $Id: lhood_site_cached_estep_shared.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __LHOOD_SITE_CACHED_ESTEP_SHARED_H
#define __LHOOD_SITE_CACHED_ESTEP_SHARED_H


#include "lhood_site_cached_shared.h"
#include "lhood_site_cached_tree_lhood_info.h"
#include "site_data_fastlup.h"
#include "subs_ml_types.h"


struct tree_walk_cache_estep_base : public tree_walk_cache_estep_top {

  tree_walk_cache_estep_base(unsigned b,
                             unsigned s)
    : tree_walk_cache_estep_top(b,s) {}

  virtual
  ~tree_walk_cache_estep_base() {}

  void
  virtual
  increment_site_stats(const tree_lhood_info& ti,
                       const site_data_fastlup& sdf,
                       const int site_id,
                       const site_data_fastlup_core::index_type* seq_state) = 0;

};


/// adaption of lhood_site_cached for em. difference is that internal
/// node info is gathered traveling both up and down each branch
///
void
lhood_site_cached_estep(const tree_lhood_info& ti,
                        const site_data_fastlup& sdf,
                        const std::vector<tree_lhood_cache>& prior_branch_lhood,
                        const tree_lhood_cache& prior_combined_lhood,
                        const int leaf_count,
                        site_data_fastlup_core::index_type* leaf_state,
                        int& site_id,
                        tree_walk_cache_estep_base& twe,
                        const bool is_reverse=false);

#endif

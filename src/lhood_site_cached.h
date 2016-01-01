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

// $Id: lhood_site_cached.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __LHOOD_SITE_CACHED_H
#define __LHOOD_SITE_CACHED_H

#include "lhood_site_cached_tree_lhood_info.h"
#include "site_data_fastlup_core.h"

#include <vector>


/// \brief subgraph-caching site lhood function (w/ site mask)
///
void
lhood_site_cached(const tree_lhood_info& ti,
                  const site_data_fastlup_core& sdf,
                  const std::vector<bool>& site_prob_mask,
                  prob_t* const site_prob);

/// \brief subgraph-caching site lhood function
///
inline
void
lhood_site_cached(const tree_lhood_info& ti,
                  const site_data_fastlup_core& sdf,
                  prob_t* const site_prob){

  std::vector<bool> local_mask(sdf.len,false);
  lhood_site_cached(ti,sdf,local_mask,site_prob);
}

#endif

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

// $Id: lhood_site_cached_estep_pprob.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __LHOOD_SITE_CACHED_ESTEP_PPROB_H
#define __LHOOD_SITE_CACHED_ESTEP_PPROB_H


#include "lhood_site_cached_tree_lhood_info.h"
#include "pprob_cache.h"
#include "site_data_fastlup.h"

/// \brief get tree node posterior probabilities (per site or all site)
///
void
lhood_site_cached_estep_pprob(const tree_lhood_info& ti,
                              const site_data_fastlup& sdf,
                              pprob_cache& ppc);

#endif

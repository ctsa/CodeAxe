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

// $Id: lhood_site_cached_optorder.h 916 2007-10-12 21:13:20Z ctsa $

/// \file

#ifndef __LHOOD_SITE_CACHED_OPTORDER_H
#define __LHOOD_SITE_CACHED_OPTORDER_H

#include "bi_tree.h"
#include "site_data_fastlup_core.h"


/// \brief chooses [nearly-]optimal leaf evaluation order for
///   subgraph caching lhood function for small trees (L<40)
///
void
lhood_site_cached_optorder(const bi_tree& tree,
                           const int n_states,
                           site_data_fastlup_core& sdf);

#endif

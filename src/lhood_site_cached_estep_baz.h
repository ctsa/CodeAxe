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

// $Id: lhood_site_cached_estep_baz.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __LHOOD_SITE_CACHED_ESTEP_BAZ_H
#define __LHOOD_SITE_CACHED_ESTEP_BAZ_H


#include "lhood_site_cached_tree_lhood_info.h"
#include "rate_gtor_nscodon_base.h"
#include "site_data_fastlup.h"

/// \brief hacky look at the bazykin-kondrashov double nonsynon mut effect
///
void
lhood_site_cached_estep_baz(const rate_gtor_nscodon_base& rg,
                            const tree_lhood_info& ti,
                            const site_data_fastlup& sdf);

#endif

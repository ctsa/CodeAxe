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

// $Id: site_data_fastlup_util.h 1173 2008-03-25 18:26:03Z ctsa $

/// \file

#ifndef __SITE_DATA_FASTLUP_UTIL_H
#define __SITE_DATA_FASTLUP_UTIL_H

#include "site_data_fastlup.h"
#include "subs_ml_model.h"

struct bi_tree;
struct site_data;


/// \brief rebuild sdf in a reduced state space
///
void
site_data_fastlup_state_reduction(const bi_tree& tree,
                                  const unsigned reduced_n_states,
                                  const site_data_fastlup& sdf_in,
                                  const unsigned* index_translation,
                                  site_data_fastlup& sdf_out);

void
mdl_site_data_fastlup_init(const site_data& sd,
                           const subs_ml_model& mdl,
                           site_data_fastlup& sdf);

#endif

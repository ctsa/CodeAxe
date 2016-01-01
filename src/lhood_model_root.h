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

// $Id: lhood_site.h 916 2007-10-12 21:13:20Z ctsa $

/// \file
/// \brief model lhood functions optimized for root-only parameter adjustment
///

#ifndef __LHOOD_MODEL_ROOT_H
#define __LHOOD_MODEL_ROOT_H

#include "simple_util.h"
#include "subs_ml_types.h"

#include <memory>

struct condition_func;
struct root_gtor;
struct site_data_fastlup;
struct subs_ml_model;


void
get_root_node_site_partial_prob_cat(const subs_ml_model& mdl,
                                    const site_data_fastlup& sdf,
                                    const unsigned cat_no,
                                    prob_t** root_spp,
                                    condition_func& cf);

void
get_root_node_site_partial_prob(const subs_ml_model& mdl,
                                const site_data_fastlup& sdf,
                                simple_matrix3d<prob_t>& root_spp,
                                std::auto_ptr<condition_func>* cf);

void
get_site_prob_from_root_spp(const root_gtor& rg,
                            const prob_t * const * root_spp,
                            const unsigned n_sites,
                            const unsigned cat_no,
                            prob_t* site_prob);

#endif

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

// $Id: cat_post_prob.h 1042 2007-12-04 01:52:18Z ctsa $

/// \file

#ifndef __CAT_POST_PROB_H
#define __CAT_POST_PROB_H


#include "cat_post_prob_mode.h"
#include "rate_gtor_options.h"
#include "site_data_fastlup.h"
#include "subs_ml_types.h"

#include <iosfwd>

struct bi_tree;
struct cat_manager;


void
make_cat_post_prob(const unsigned n_cats,
                   const site_data_fastlup& full_sdf,
                   const CPPM::index_t cat_post_prob_mode,
                   prob_t** ppi);

void
cat_post_prob_report(const RATE_GTOR_MODEL::index_t rgm,
                     const bi_tree& tree,
                     const cat_manager& cm,
                     const site_data_fastlup& full_sdf,
                     const prob_t * const * ppi,
                     const CPPM::index_t cat_post_prob_mode,
                     std::ostream& os);

#endif

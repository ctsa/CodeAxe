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

// $Id: lhood_model_em_baz.h 1149 2008-03-13 00:46:06Z ctsa $

/// \file

#ifndef __LHOOD_MODEL_EM_BAZ_H
#define __LHOOD_MODEL_EM_BAZ_H


struct bi_tree;
struct site_data_fastlup;
struct subs_ml_model;
struct tree_probs;

/// estep update
///
/// !!ALL TPROBS TRANSPOSED!!
void
update_etrans_baz(const subs_ml_model& mdl,
                  const bi_tree& tree,
                  const tree_probs& tprob,
                  const site_data_fastlup& sdf,
                  const unsigned n_states);

#endif

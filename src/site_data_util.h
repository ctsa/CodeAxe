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

// $Id: site_data_util.h 1163 2008-03-21 00:13:32Z ctsa $

/// \file

#ifndef __SITE_DATA_UTIL_H
#define __SITE_DATA_UTIL_H

struct site_data;
struct nuc_seq_data;

#include "rate_gtor_options.h"

/// \brief convert continuous nucleotide data to an unordered set of
/// sites
///
/// \param is_no_adjacent_nuc_diff filter sites at which adjacent
///                                nucleotides columns are not
///                                constant across all taxa
///
/// \param is_single_nuc_diff filter sites at which more than one
///                           nucleotide column is not constant across
///                           all taxa
///
void
nuc_seq_to_site_data(site_data& sd,
                     const nuc_seq_data& in_seq,
                     const RATE_GTOR_MODEL::index_t model,
                     const bool is_codon_border,
                     const bool is_nuc_4fold_filter,
                     const double gc_max = 1.,
                     const double gc_min = -1.,
                     const bool is_no_adjacent_nuc_diff = false,
                     const bool is_single_nuc_diff = false);
#endif

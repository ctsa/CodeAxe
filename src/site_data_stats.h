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

// $Id: site_data_stats.h 1170 2008-03-24 20:55:56Z ctsa $

/// \file

#ifndef __SITE_DATA_STATS_H
#define __SITE_DATA_STATS_H

#include "subs_ml_types.h"

struct site_data;

#include <vector>


/// \brief get symbol distribution from site data for all taxa and
/// data classes
///
/// \param distro [n_states]
///
void
site_data_2_state_pdf(const site_data& sd,
                      const unsigned n_states,
                      prob_t* distro); // [n_states]

void
site_data_2_taxa_state_count(const site_data& sd,
                             const unsigned n_states,
                             unsigned* count, // [n_orgs*n_states]
                             const bool is_accumulate = false,
                             const bool is_single_data_class = false,
                             const unsigned target_data_class_id = 0);

void
site_data_2_taxa_state_pdf(const site_data& sd,
                           const unsigned n_states,
                           prob_t* distro);  // [n_orgs*n_states]

void
site_data_2_data_class_state_count(const site_data& sd,
                                   const unsigned n_data_classes,
                                   const unsigned n_states,
                                   std::vector<unsigned>& class_state_count);  // [n_data_classes*n_states]

void
site_data_2_data_class_state_pdf(const site_data& sd,
                                 const unsigned n_data_classes,
                                 const unsigned n_states,
                                 prob_t* distro);  // [n_data_classes*n_states]

void
site_data_2_data_class_count(const site_data& sd,
                             const unsigned n_data_classes,
                             std::vector<unsigned>& class_count);  // [n_data_classes]

void
site_data_counts(const site_data& sd,
                 unsigned& site_count,
                 unsigned& conserved_site_count,
                 unsigned& unique_site_count);

#endif

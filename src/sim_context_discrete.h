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

// $Id: sim_context_discrete.h 1101 2008-01-21 23:44:08Z ctsa $

/// \file

#ifndef __SIM_CONTEXT_DISCRETE_H
#define __SIM_CONTEXT_DISCRETE_H


#include "subs_ml_types.h"

#include <vector>


/// \brief discrete sequence simulation along a single tree branch
///
/// \param startseq input sequence
/// \param endeq output sequence
/// \param time branch time
/// \param unit_time time segment used for discrete evolutionary steps
///
void
simulate_discrete_time_branch_context(const std::vector<unsigned>& startseq,
                                      const std::vector<unsigned>& cat_seq,
                                      const std::vector<unsigned>& group_seq,
                                      std::vector<unsigned>& endseq,
                                      const smlfloat time,
                                      const smlfloat unit_time,
                                      const prob_t* const * sitesub_prob,
                                      const prob_t* const * sitesub_cdf,
                                      const prob_t * const * bg_nuc_cdf_cat,
                                      std::vector<smlfloat>& cat_sim_time,
                                      const bool is_codon_site_model,
                                      const bool is_report_sim_time,
                                      const prob_t* const * sitesub_prob_neutral);

#endif

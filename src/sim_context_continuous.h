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

// $Id: sim_context_continuous.h 1130 2008-01-29 03:41:42Z ctsa $

/// \file

#ifndef __SIM_CONTEXT_CONTINUOUS_H
#define __SIM_CONTEXT_CONTINUOUS_H


#include "subs_ml_types.h"

#include <vector>


///
void
simulate_continuous_time_branch_context(const std::vector<unsigned>& startseq,
                                        const std::vector<unsigned>& cat_seq,
                                        const std::vector<unsigned>& group_seq,
                                        std::vector<unsigned>& endseq,
                                        const smlfloat time,
                                        const unsigned n_cats,
                                        const prob_t* const * sitesub_prob,
                                        const prob_t* const * sitesub_cdf,
                                        const prob_t* const * bg_nuc_cdf_cat,
                                        std::vector<smlfloat>& cat_sim_time,
                                        const bool is_codon_site_model,
                                        const bool is_report_sim_time,
                                        const prob_t* const * sitesub_prob_neutral);

#endif

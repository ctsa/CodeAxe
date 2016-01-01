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

// $Id: sim_util.h 1193 2008-03-29 03:14:54Z ctsa $

/// \file

#ifndef __SIM_UTIL_H
#define __SIM_UTIL_H

#include "nuc_seq_data.h"
#include "subs_ml_model.h"
#include "site_data.h"

#include <iosfwd>
#include <vector>

struct cat_manager;
struct sim_options;
struct subs_ml_model;


extern const smlfloat SIMULATE_DISCRETE_TIME_UNIT;


void
get_ancestral_seq(const subs_ml_model& mdl,
                  const unsigned size,
                  const std::vector<unsigned>& cat_seq,
                  std::vector<unsigned>& seq);

void
get_sim_group_seq(std::vector<unsigned>& group_seq,
                  const sim_options& sim_opt);

void
get_sim_cat_seq(std::vector<unsigned>& cat_seq,
                const std::vector<unsigned>& group_seq,
                const unsigned size,
                const subs_ml_model& mdl);

#if 0
void
sim_site_data_init(const bi_tree& tree,
                   const RATE_GTOR_MODEL::index_t rgm,
                   site_data& sd);
#endif

void
sim_seq_data_init(const bi_tree& tree,
                  const unsigned n_groups,
                  const unsigned n_cats,
                  const char * const gtag,
                  const cat_manager& cm,
                  nuc_seq_data& nsd);

void
report_sim_tree_times(const cat_manager& cm,
                      const bi_tree& tree,
                      const std::vector<std::vector<smlfloat> >& cat_sim_time,
                      std::ostream& os);


struct group_range {

  group_range(const std::vector<unsigned>& group_seq,
              const unsigned left_pad=0,
              const unsigned right_pad=0)
    : _gs(group_seq), _lp(left_pad), _rp(right_pad), _lgs(0) {}

  /// must call in group order && group_seq must be sorted
  bool
  get_gbounds(const unsigned group_no,
              unsigned& gstart,
              unsigned& gstop);

private:
  const std::vector<unsigned>& _gs;
  const unsigned _lp;
  const unsigned _rp;
  unsigned _lgs;
};


#endif

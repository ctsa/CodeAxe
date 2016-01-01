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

// $Id: sim_options.h 1070 2007-12-24 21:12:35Z ctsa $

/// \file

#ifndef __SIM_OPTIONS_H
#define __SIM_OPTIONS_H

#include "rate_gtor_options.h"
#include "subs_ml_types.h"


namespace SIM_MODEL {
  enum index_t { CONTINUOUS,
                 DISCRETE,
                 ISS };
}


struct sim_options {
  sim_options() :
    size(0),
    method(SIM_MODEL::CONTINUOUS),
    is_report_time(false),
    group_size_min(100),
    group_size_max(200),
    assigned_cat_prob(0.),
    output_site_model(RATE_GTOR_MODEL::NONE),
    is_nuc_seq_output(false) {}

  unsigned size;
  SIM_MODEL::index_t method;
  bool is_report_time;
  unsigned group_size_min;
  unsigned group_size_max;
  prob_t assigned_cat_prob;
  RATE_GTOR_MODEL::index_t output_site_model;
  bool is_nuc_seq_output;
};


#endif

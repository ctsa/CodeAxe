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

// $Id: subs_ml_confidence.h 1065 2007-12-14 21:49:28Z ctsa $

/// \file

#ifndef __SUBS_ML_CONFIDENCE_H
#define __SUBS_ML_CONFIDENCE_H


#include "site_data_fastlup.h"
#include "subs_ml_model.h"

#include <string>
#include <vector>


const smlfloat DEFAULT_ALPHA(0.05);


struct conf_options {

  conf_options()
    : alpha(DEFAULT_ALPHA), is_fisher_method(false), is_single_param(false) {}

  double alpha;
  bool is_fisher_method;
  bool is_single_param;
};


void
get_subs_ml_model_conf_interval(const site_data_fastlup& sdf,
                                const subs_ml_model& mdl,
                                const conf_options& co,
                                const simple_array<bool>& ci_mask,
                                const std::string& outfile);

#endif

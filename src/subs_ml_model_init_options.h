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

// $Id: subs_ml_model_init_options.h 1117 2008-01-28 19:48:39Z ctsa $

/// \file

#ifndef __SUBS_ML_MODEL_INIT_OPTIONS_H
#define __SUBS_ML_MODEL_INIT_OPTIONS_H


#include "bg_gtor_options.h"
#include "cat_manager.h"
#include "obs_info.h"
#include "param_init_type.h"
#include "rate_gtor_options.h"
#include "rate_gtor_nuc_options.h"
#include "rate_gtor_nscodon_options.h"
#include "root_gtor_options.h"
#include "time_gtor.h"

#include "util/general/uncopyable.h"

#include <string>

struct audit_info;



struct subs_ml_model_init_options : private uncopyable {

  subs_ml_model_init_options(const RATE_GTOR_MODEL::index_t init_ragm,
                             const audit_info& init_ai)
    : ragm(init_ragm),
      c5t(C5_APPROX::NONE),
      rogm(ROOT_GTOR_MODEL::FULL),
      pinit(PARAM_INIT_TYPE::START),
      ai(init_ai) {}


  std::string tree_buffer;
  time_gtor_options tgo;
  const RATE_GTOR_MODEL::index_t ragm;
  rate_gtor_options ropt;
  rate_gtor_nuc_options nopt;
  rate_gtor_nscodon_options copt;
  bg_gtor_options bopt;
  C5_APPROX::index_t c5t;
  ROOT_GTOR_MODEL::index_t rogm;
  cat_manager_options catman_opt;
  PARAM_INIT_TYPE::index_t pinit;
  const audit_info& ai;
};

#endif

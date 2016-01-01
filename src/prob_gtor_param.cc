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

// $Id: prob_gtor_param.cc 1113 2008-01-27 00:23:59Z ctsa $

/// \file

#include "prob_gtor_param.h"
#include "subs_ml_print_util.h"

#include <iomanip>
#include <iostream>


void
prob_gtor_param::
store_state(const char* label,
            std::ostream& os) const {

  os << std::setprecision(SUBS_ML_PRINT_MODEL_PRECISION);

  const unsigned ps(param_size());
  for(unsigned j(0);j<ps;++j){
    os << label << "_param_" << j << " " << _is_train_param[j] << " " << _param[j] << "\n";
  }
}

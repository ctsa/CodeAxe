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

// $Id: rate_gtor_util.cc 885 2007-10-02 02:40:44Z ctsa $

/// \file

#include "rate_gtor_util.h"
#include "stationary_pdistro.h"
#include "subs_ml_ptol.h"
#include "util/general/log.h"
#include "util/math/prob_util_io.h"

#include <cstdlib>

#include <exception>
#include <ostream>


void
get_stationary_pdistro(smlfloat* stat_pdistro,
                       const rate_gtor& rg,
                       const rates_func_options& opt) {

  const unsigned ss(rg.state_size());
  smlfloat* rn(new smlfloat[ss*ss]);
  rg.rates(rn,opt);
  try{
    get_stationary_pdistro_from_rates(stat_pdistro,rn,ss);
  } catch (std::exception& e){
    log_os << "EXCEPTION: " << e.what() << "\n"
           << "...caught in " << __FILE__ << ":" << __LINE__ << ": dumping rate_gtor state: \n";
    rg.store_state(log_os);
    throw;
  }

  delete [] rn;
#ifdef DEBUG
  pdistro_check(stat_pdistro,ss,SUBS_ML_PTOL);
#endif
}

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

// $Id: prob_util.cc 871 2007-10-01 00:17:23Z ctsa $

/// \file

#include "prob_util.h"
#include "../general/log.h"

#include <cstdlib>

#include <ostream>
#include <iomanip>


void
pdistro_check_fail_val(const double v,
                       const unsigned i){
  log_os << "pdistro_check: Negative/NAN/inf p in prob distro: p[i],i "
         << v << " , " << i << "\n";
  abort();
}




void
pdistro_check_fail_sum(const double s){

  log_os << "pdistro_check: bad distro sum: "
         << std::setprecision(12) << s << "\n";
  abort();
}

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

// $Id: subs_ml_print_util.cc 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "subs_ml_print_util.h"

#include <ostream>
#include <iomanip>

void
lnp_line_report(const smlfloat lnp,
                const smlfloat lnp_norm,
                const unsigned site_count,
                std::ostream& os){
  os << std::setprecision(SUBS_ML_PRINT_LNP_PRECISION)
     << std::setw(14) << lnp << " "
     << std::setw(14) << lnp_norm << " "
     << std::setw(14) << lnp/site_count << "\n";
}

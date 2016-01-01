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

// $Id: indy_random_var_io.hh 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "chi_sqr.h"

#include <cmath>

#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>


template <typename T>
std::ostream&
operator<<(std::ostream& os,
           const irv_t<T>& p) {

  static const T alpha(0.05);
  static const T lnp_conf_diff(static_cast<T>(chi_sqr(1.-alpha,1.)/2.));

  if(p.v != 0.){
    unsigned pr(os.precision());
    std::ostringstream oss;
    oss.precision(pr);
    oss << p.m;
    oss << " +/- ";
    const T target_delta_plus(std::sqrt(p.v*(2.*lnp_conf_diff)));
    oss << target_delta_plus;
    os << oss.str();
  } else {
    os << p.m;
  }

  return os;
}


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

// $Id: report_util.h 916 2007-10-12 21:13:20Z ctsa $

/// \file

#ifndef __REPORT_UTIL_H
#define __REPORT_UTIL_H


#include "bi_tree.h"
#include "subs_ml_types.h"
#include "util/math/indy_random_var.h"


#include <iosfwd>
#include <vector>

template <typename FloatType>
void
report_time_instance(const char* time_label,
                     const std::vector<FloatType>& branch_time,
                     const bi_tree& tree,
                     std::ostream& os);

#include "report_util.hh"

#endif

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

// $Id: bioseq_report.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __BIOSEQ_REPORT_H
#define __BIOSEQ_REPORT_H


#include "subs_ml_types.h"

#include <iosfwd>


template <typename FloatType>
void
report_sorted_nsaa_pair_values(const FloatType* nsaa_val,
                               std::ostream& os,
                               const char* label);

#include "bioseq_report.hh"

#endif

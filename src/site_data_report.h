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

// $Id: site_data_report.h 1086 2008-01-11 20:23:09Z ctsa $

/// \file

#ifndef __SITE_DATA_REPORT_H
#define __SITE_DATA_REPORT_H

#include "subs_ml_types.h"

#include <iosfwd>

struct site_data;
struct subs_ml_model;

void
site_data_report(const site_data& sd,
                 const bool is_pretty_print_data,
                 std::ostream& os);

/// temporary!!
void
site_data_cat_split(const site_data& sd,
                    const char* outtag);

#if 0
void
report_nscodon_model(const subs_ml_model& mdl,
                     std::ostream& os);

void
report_multinuc_model(const subs_ml_model& mdl,
                   std::ostream& os);

void
report_trinuc_model(const subs_ml_model& mdl,
                    std::ostream& os);
#endif

#endif

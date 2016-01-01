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

// $Id: rate_gtor_nuc_base_util.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __RATE_GTOR_NUC_BASE_UTIL_H
#define __RATE_GTOR_NUC_BASE_UTIL_H


#include "rate_gtor_nuc_base.h"


void
rates_nuc_context(const rate_gtor_nuc_base& r,
                  smlfloat rates[],
                  const rates_func_options& opt);

#endif

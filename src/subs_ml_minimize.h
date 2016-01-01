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

// $Id: subs_ml_minimize.h 1179 2008-03-26 00:23:59Z ctsa $

/// \file

#ifndef __SUBS_ML_MINIMIZE_H
#define __SUBS_ML_MINIMIZE_H


#include "estep_dat_t.h"
#include "site_data_fastlup.h"
#include "subs_ml_model.h"


/// "expert" version
///
void
sml_fullstep_min_tol(const site_data_fastlup& sdf,
                     subs_ml_model& mdl,
                     const smlfloat start_ntol,
                     const smlfloat end_ntol,
                     const bool is_skip_hybrid_cg = false);


void
sml_fullstep_min(const site_data_fastlup& sdf,
                 subs_ml_model& mdl,
                 const smlfloat start_lnp,
                 const bool is_start_tol = true);

#if 0
void
sml_mstep_min(const estep_dat_t& edat,
              subs_ml_model& mdl);
#endif

void
sml_root_cycle_min(const site_data_fastlup& sdf,
                   subs_ml_model& mdl,
                   const smlfloat start_lnp);

#endif

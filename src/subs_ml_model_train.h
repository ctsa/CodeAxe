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

// $Id: subs_ml_model_train.h 1022 2007-11-13 22:09:39Z ctsa $

/// \file

#ifndef __SUBS_ML_MODEL_TRAIN_H
#define __SUBS_ML_MODEL_TRAIN_H

#include "site_data_fastlup.h"
#include "subs_ml_model.h"

#if 0
void
train_subs_ml_model_obs_partition(const site_data_fastlup& sdf,
                                  subs_ml_model& mdl);
#endif

/// \brief get parameter mle's for site data
///
void
train_subs_ml_model(const site_data_fastlup& sdf,
                    subs_ml_model& mdl,
                    bool is_refine = false); // tmp shoe-in control

#endif

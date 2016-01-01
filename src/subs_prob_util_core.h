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

// $Id: subs_prob_util_core.h 1055 2007-12-08 04:27:31Z ctsa $

/// \file

#ifndef __SUBS_PROB_UTIL_CORE_H
#define __SUBS_PROB_UTIL_CORE_H

#include "subs_ml_model.h"
#include "util/general/workspace.h"


/// \brief the 'expert' (ugly) version of subs_prob_util
///
/// calculates subs prob for a single opt and a series of times, when
/// that time position is not masked
///
/// accepts any iterator of times assumed to be over the same rate matrix,
/// (ie. 'opt' applies to all branches)
///
template <typename RandomAccessIterator>
void
get_subs_prob(prob_t** subs_prob,
              const subs_ml_model& mdl,
              const rates_func_options& opt,
              const unsigned n_subs_prob,
              RandomAccessIterator subs_prob_time,
              workspace<char> ws[3],
              const bool* is_mask=0);

#include "subs_prob_util_core.hh"

#endif

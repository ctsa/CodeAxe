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

// $Id: cat_type.cc 1153 2008-03-18 00:24:07Z ctsa $

/// \file

#include "cat_type.h"
#include "util/general/die.h"

#include <cassert>


const char* CAT_PARAM_TYPE::syms[] =
  {"mut_rate","mut_model","sel_strength","sel_matrix","root","time","obs","subs_rate_bg"};
const char* CAT_PARAM_TYPE::short_syms[] =
  {"mr","mm","ss","sm","ro","ti","ob","bg"};

const bool CAT_PARAM_TYPE::is_mixable[] =
  { true, false, true, false, false, false, false, false };

const char* CAT_MIX_TYPE::syms[] = {"site","group"};
const char* CAT_MIX_TYPE::short_syms[] = {"s","g"};


//const unsigned CAT_TYPE_SIZE(CAT_PARAM_TYPE::SIZE*CAT_MIX_TYPE::SIZE);

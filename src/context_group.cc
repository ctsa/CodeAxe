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

// $Id: context_group.cc 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "context_group.h"
#include "util/general/die.h"


CONTEXT_GROUP::index_t
context_model_group(const CONTEXT_MODEL_NUC::index_t cmn){

  using namespace CONTEXT_MODEL_NUC;

  switch(cmn){
  case INDY:         return CONTEXT_GROUP::NONE;
  case PRE_DOUBLET:  return CONTEXT_GROUP::ONLY_PRE;
  case POST_DOUBLET: return CONTEXT_GROUP::ONLY_POST;
  case FACTORED_TRIPLET:
  case TRIPLET:
  case CPG_ONLY:
  case TPA_ONLY:
  case CPG_1TI:
  case CPG_2TI:
  case CPG_NONREV:   return CONTEXT_GROUP::TWOWAY;
  default :
    die("Unsupported nuc context model");
    return CONTEXT_GROUP::NONE;
  }
}

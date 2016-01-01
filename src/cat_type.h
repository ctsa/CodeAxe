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

// $Id: cat_type.h 1107 2008-01-25 02:40:20Z ctsa $

/// \file

#ifndef __CAT_TYPE_H
#define __CAT_TYPE_H

#include <string>


namespace CAT_PARAM_TYPE {
  enum index_t {
    MUT_RATE,
    MUT_MODEL,
    SEL_STRENGTH,
    SEL_MATRIX,
    ROOT,
    TIME,
    OBS,
    SUBS_RATE_BG,
    SIZE};

  extern const char* syms[SIZE];
  extern const char* short_syms[SIZE];

  extern const bool is_mixable[SIZE];
}

// this enum cannot be extended
//
namespace CAT_MIX_TYPE {
  enum index_t {
    EITHER=-1,
    SITE,
    GROUP,
    SIZE
  };

  extern const char* syms[SIZE];
  extern const char* short_syms[SIZE];
}


enum { CAT_TYPE_SIZE=CAT_PARAM_TYPE::SIZE*CAT_MIX_TYPE::SIZE };



inline
unsigned
typed_cat_index(const CAT_PARAM_TYPE::index_t pt,
                const CAT_MIX_TYPE::index_t mt) {
  return mt*CAT_PARAM_TYPE::SIZE+pt;
}


#endif

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

// $Id: site_data_fastlup_core_util.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __SITE_DATA_FASTLUP_CORE_UTIL_H
#define __SITE_DATA_FASTLUP_CORE_UTIL_H

#include "site_data_fastlup_core.h"

/// map the order of entries in sdf to the location of their
/// corresponding index in sdf_superset. if is_translate, then remap
/// sdf's index states before searching sdf_superset for the
/// corresponding remapped index
///
template <typename RandomAccessIterator,
          typename OutputIterator>
void
map_translated_index_superset_data(const site_data_fastlup_core& sdf,
                                   const site_data_fastlup_core& sdf_superset,
                                   RandomAccessIterator index_translation,
                                   OutputIterator index_map,
                                   const bool is_translate=true);

/// map the order of entries in sdf to the location of their
/// corresponding index in sdf_superset.
///
template <typename OutputIterator>
void
map_superset_data(const site_data_fastlup_core& sdf,
                  const site_data_fastlup_core& sdf_superset,
                  OutputIterator index_map){
  map_translated_index_superset_data(sdf,sdf_superset,0,index_map,false);
}

#include "site_data_fastlup_core_util.hh"

#endif

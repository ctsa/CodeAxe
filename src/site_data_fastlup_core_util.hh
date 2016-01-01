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

// $Id: site_data_fastlup_core_util.hh 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "util/general/die.h"

#include <algorithm>


template <typename RandomAccessIterator,
          typename OutputIterator>
void
map_translated_index_superset_data(const site_data_fastlup_core& sdf,
                                   const site_data_fastlup_core& sdf_superset,
                                   RandomAccessIterator index_translation,
                                   OutputIterator index_map,
                                   const bool is_translate){

  typedef const site_data_fastlup_cell_core* iter_t;

  if(sdf.n_orgs != sdf_superset.n_orgs){ die("conflicting sdf org counts\n"); }

  iter_t superset_begin(sdf_superset.data_core);
  iter_t superset_end(sdf_superset.data_core+sdf_superset.len);

  site_data_fastlup_cell_core search_cell;
  if(is_translate){
    search_cell.index=new site_data_fastlup_cell_core::index_type[sdf.n_orgs];
  }

  const site_data_fastlup_cell_core* search_cell_ptr(&search_cell);

  for(unsigned i(0);i<sdf.len;++i){
    if(is_translate){
      for(unsigned t(0);t<sdf.n_orgs;++t){
        search_cell.index[t] = *(index_translation+sdf.data_core[i].index[t]);
      }
    } else {
      search_cell_ptr = &(sdf.data_core[i]);
    }

    iter_t j(std::lower_bound(superset_begin,superset_end,
                              *search_cell_ptr,sdf_superset.cell_sorter));
    if(j==superset_end){ die("invalid sdf_superset in map_superset_data"); }
    *(index_map++) = j-superset_begin;
  }

  if(is_translate) delete [] search_cell.index;
}

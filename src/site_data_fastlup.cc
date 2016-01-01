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

// $Id: site_data_fastlup.cc 1182 2008-03-27 01:55:02Z ctsa $

/// \file

#include "site_data_fastlup.h"
#include "util/general/die.h"

#include <algorithm>



void
site_data_fastlup::
init(unsigned init_len,
     unsigned init_n_orgs,
     unsigned init_n_assigned_data_sets) {

  clear(false);
  base_t::init(init_len,init_n_orgs);

  data = new site_data_fastlup_cell[len];

  n_assigned_data_sets = init_n_assigned_data_sets;

  _ads_site_mask.resize(n_assigned_data_sets);
  for(unsigned i(0);i<n_assigned_data_sets;++i){ _ads_site_mask[i].resize(len,true); }

  for(unsigned i(0);i<len;++i){
    data[i].index=cell_index_pool+i*n_orgs;
  }
}



void
site_data_fastlup::
clear(bool is_clear_base) {
  _ads_site_mask.clear();
  if(data) delete [] data; data=0;
  if(is_clear_base) base_t::clear();
}



void
site_data_fastlup::
fix(){
  typedef site_data_fastlup_cell::group_assigned_data_set_vector::const_iterator gavi;

  if(_is_fixed) { die("second fix in site_data_fastlup"); }

  /// \todo gcc is using more memory than expected on this sort call,
  /// apparently because there's data in cell_sorter
  ///
  std::sort(data,data+len,cell_sorter);

  // Calculate ads_mask, n_groups and total_count. Calculate n_groups
  // assuming group numbers have already been compressed:
  n_groups=0;
  total_count=0;
  for(unsigned i(0);i<len;++i){
    gavi g(data[i].group_assigned_data_set_count.begin());
    const gavi g_end(data[i].group_assigned_data_set_count.end());
    for(;g!=g_end;++g){
      const site_data_fastlup_code::group_id_type group_id(g->first.group_id);
      const site_data_fastlup_code::assigned_data_set_id_type adset_id(g->first.assigned_data_set_id);
      const unsigned count(g->second);
      n_groups = std::max(n_groups,static_cast<unsigned>(group_id+1));
      total_count += count;
      _ads_site_mask[adset_id][i] = false;
    }
  }

  _is_fixed=true;

  base_t::fix();
}

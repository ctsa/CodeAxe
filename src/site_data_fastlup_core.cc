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

// $Id: site_data_fastlup_core.cc 1182 2008-03-27 01:55:02Z ctsa $

/// \file

#include "site_data_fastlup_core.h"
#include "util/general/die.h"

#include <algorithm>


site_data_fastlup_core&
site_data_fastlup_core::
operator=(const site_data_fastlup_core& rhs){
  if(&rhs == this) return *this;
  clear();
  len = rhs.len;
  n_orgs = rhs.n_orgs;
  data_core = new site_data_fastlup_cell_core[len];
  is_block_boundary = rhs.is_block_boundary;
  order = rhs.order;
  cell_index_pool = new index_type[len*n_orgs];
  for(unsigned i(0);i<(len*n_orgs);++i){
    cell_index_pool[i] = rhs.cell_index_pool[i];
  }

  for(unsigned i(0);i<len;++i){
    const unsigned offset(rhs.data_core[i].index-rhs.cell_index_pool);
    data_core[i].index=cell_index_pool+offset;
  }

  return *this;
}


void
site_data_fastlup_core::
init(unsigned init_len,
     unsigned init_n_orgs) {

  clear();
  len = init_len;
  if(len==0) { die("Invalid len in sdf init"); }

  data_core = new site_data_fastlup_cell_core[len];
  is_block_boundary.resize(len);

  n_orgs = init_n_orgs;
  if(n_orgs==0) { die("Invalid n_org in sdf init"); }

  order.resize(n_orgs);
  for(unsigned i(0);i<n_orgs;++i){ order[i] = i; }

  cell_index_pool = new index_type[len*n_orgs];
  for(unsigned i(0);i<len;++i){
    data_core[i].index=cell_index_pool+i*n_orgs;
  }
}




void
site_data_fastlup_core::
clear() {

  if(len) {
    if(cell_index_pool) delete [] cell_index_pool;
    cell_index_pool = 0;
    n_orgs = 0;

    delete [] data_core;
    is_block_boundary.resize(0);
    len = 0;
  }
}




void
site_data_fastlup_core::
fix(){

  if(_is_fixed) die("double fixing sdf_core");

  std::sort(data_core,data_core+len,cell_sorter);

  const unsigned last_org(order.back());
  std::vector<unsigned> last_index(n_orgs);

  for(unsigned i(0);i<len;++i){
    is_block_boundary[i]=false;
    for(unsigned j(0);j<n_orgs;++j){
      if(j==last_org) continue;
      if(i==0 || data_core[i].index[j] != last_index[j]){
        is_block_boundary[i]=true;
        last_index[j] = data_core[i].index[j];
      }
    }
  }

  _is_fixed=true;
}

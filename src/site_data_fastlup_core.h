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

// $Id: site_data_fastlup_core.h 1182 2008-03-27 01:55:02Z ctsa $

/// \file

#ifndef __SITE_DATA_FASTLUP_CORE_H
#define __SITE_DATA_FASTLUP_CORE_H

#include "site_code.h"

#include <vector>


struct site_data_fastlup_cell_core {

  typedef site_code::index_type index_type;

  index_type* index;
};


struct cell_compare_func {

  explicit
  cell_compare_func(const std::vector<unsigned>& f) : order(f) {}

  bool operator()(const site_data_fastlup_cell_core& lhs,
                  const site_data_fastlup_cell_core& rhs) const {

    const unsigned n_index(order.size());
    for(unsigned i(0);i<n_index;++i){
      const unsigned oi(order[i]);
      if(lhs.index[oi] != rhs.index[oi]){
        return lhs.index[oi] < rhs.index[oi];
      }
    }
    return false;
  }

  const std::vector<unsigned>& order;
};




/// \brief a fast access site data structure
///
/// this object is only a listing of site indices -- the count, group
/// and assigned data set info are tracked in derived class
///
struct site_data_fastlup_core {
private:
  site_data_fastlup_core(const site_data_fastlup_core&);

public:
  typedef site_data_fastlup_cell_core::index_type index_type;

  site_data_fastlup_core()
    : len(0),n_orgs(0),
      data_core(0),order(),cell_sorter(order),cell_index_pool(0), _is_fixed(false) {}

  ~site_data_fastlup_core() { clear(); };

  site_data_fastlup_core& operator=(const site_data_fastlup_core& rhs);

  void
  init(unsigned init_len,
       unsigned init_n_orgs);

  /// sort and setup for fast indexing
  void fix();

protected:
  void clear();

public:
  unsigned len;
  unsigned n_orgs;

  site_data_fastlup_cell_core* data_core;
  std::vector<bool> is_block_boundary;

  std::vector<unsigned> order;

  cell_compare_func cell_sorter;

protected:
  index_type* cell_index_pool;
private:
  bool _is_fixed;
};


#endif

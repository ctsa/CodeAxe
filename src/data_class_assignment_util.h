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

// $Id: data_class_assignment_util.h 1222 2008-05-22 23:10:06Z ctsa $

/// \file
///
/// the cat_manager holds a mapping from data class labels to assigned
/// data sets, to which parameterized model categories are tied. this
/// structure creaties a mapping from data class ids to assigned data
/// set ids for a specific instance of site_data
///

#ifndef __DATA_CLASS_ASSIGNMENT_UTIL_H
#define __DATA_CLASS_ASSIGNMENT_UTIL_H

#include "substk_exception.h"

#include <vector>


struct site_data;
struct cat_manager;



struct data_class_id_assignment_map {

  data_class_id_assignment_map(const cat_manager& cm,
                               const site_data& sd);

  /// \todo fix site_data_fastlup state reducer so that this 'naive'
  /// ctor is not required:
  data_class_id_assignment_map(const unsigned n)
    : _is_dci_ads_map(n,true), _dci_ads_map(n), _assigned_data_set_size(n) {
    for(unsigned i(0);i<n;++i) _dci_ads_map[i] = i;
  }

  bool
  is_data_class_assigned(const unsigned data_class_id) const {
    if(data_class_id >= _dci_ads_map.size()){
      throw substk_exception("data_class_id_assignment_map.is_data_class_assigned(): invalid data_class_id");
    }

    return _is_dci_ads_map[data_class_id];
  }

  unsigned
  get_data_class_assigned_data_set_id(const unsigned data_class_id) const {
    if(! is_data_class_assigned(data_class_id)){
      throw substk_exception("data_class_id_assignment_map.get_data_class_assigned_data_set_id(): no assigned data set for data class id");
    }

    return _dci_ads_map[data_class_id];
  }

  unsigned
  assigned_data_set_size() const { return _assigned_data_set_size; }

private:
  std::vector<bool> _is_dci_ads_map;
  std::vector<unsigned> _dci_ads_map;
  unsigned _assigned_data_set_size;
};


#endif

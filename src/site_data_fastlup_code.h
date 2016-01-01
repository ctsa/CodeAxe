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

// $Id: site_data_fastlup_code.h 1165 2008-03-21 03:32:14Z ctsa $

/// \file

#ifndef __SITE_DATA_FASTLUP_CODE_H
#define __SITE_DATA_FASTLUP_CODE_H


#include "util/general/die.h"

#include <limits>


/// \brief stores all info (group_id, assigned_data_set_id) associated
/// with fastlup count information at each site and sorts
/// appropriately for c++ stdlib use
///
struct site_data_fastlup_code {

  typedef unsigned short group_id_type;
  typedef unsigned char assigned_data_set_id_type;

  explicit
  site_data_fastlup_code(unsigned _gi = 0,
                         unsigned _ai = 0)
    : group_id(_gi), assigned_data_set_id(_ai) {
    if(_gi > std::numeric_limits<group_id_type>::max())
      die("count_data_code(): group_id exceeds bounds of numerical type.");
    if(_ai > std::numeric_limits<assigned_data_set_id_type>::max())
      die("count_data_code(): assigned_data_set_id exceeds bounds of numerical type");
  }

  bool
  operator==(const site_data_fastlup_code& right) const {
    if(! (group_id == right.group_id)) return false;
    if(! (assigned_data_set_id == right.assigned_data_set_id)) return false;
    return true;
  }

  bool
  operator<(const site_data_fastlup_code& right) const {
    if(! (group_id == right.group_id)) return group_id < right.group_id;
    if(! (assigned_data_set_id == right.assigned_data_set_id)) return assigned_data_set_id < right.assigned_data_set_id;
    return false;
  }

  group_id_type group_id;
  assigned_data_set_id_type assigned_data_set_id;
};


#endif

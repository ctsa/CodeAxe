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

// $Id: site_data_code.h 1164 2008-03-21 03:28:49Z ctsa $

/// \file

#ifndef __SITE_DATA_CODE_H
#define __SITE_DATA_CODE_H


#include "util/general/die.h"

#include <limits>


/// \brief stores all info (group_id, data_class_id) associated with input
/// count information at each site and sorts appropriately for c++
/// stdlib use
///
struct site_data_code {

  typedef unsigned short group_id_type;
  typedef unsigned char data_class_id_type;

  explicit
  site_data_code(unsigned _gi = 0,
                 unsigned _di = 0)
    : group_id(_gi), data_class_id(_di) {
    if(_gi > std::numeric_limits<group_id_type>::max())
      die("site_data_code(): group_id exceeds bounds of numerical type.");
    if(_di > std::numeric_limits<data_class_id_type>::max())
      die("site_data_code(): data_class_id exceeds bounds of numerical type");
  }

  bool
  operator==(const site_data_code& right) const {
    if(! (group_id == right.group_id)) return false;
    if(! (data_class_id == right.data_class_id)) return false;
    return true;
  }

  bool
  operator<(const site_data_code& right) const {
    if(! (group_id == right.group_id)) return group_id < right.group_id;
    if(! (data_class_id == right.data_class_id)) return data_class_id < right.data_class_id;
    return false;
  }

  group_id_type group_id;
  data_class_id_type data_class_id;
};


#endif

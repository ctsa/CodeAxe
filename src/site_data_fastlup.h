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

// $Id: site_data_fastlup.h 1182 2008-03-27 01:55:02Z ctsa $

/// \file

#ifndef __SITE_DATA_FASTLUP_H
#define __SITE_DATA_FASTLUP_H

#include "site_data_fastlup_code.h"
#include "site_data_fastlup_core.h"
#include "subs_ml_types.h"
#include "util/general/uncopyable.h"

#include <string>
#include <vector>
#include <utility>


struct site_data_fastlup_cell : public site_data_fastlup_cell_core {

  typedef std::vector<std::pair<site_data_fastlup_code,unsigned> > group_assigned_data_set_vector; // (group,adset), count

  group_assigned_data_set_vector group_assigned_data_set_count;
};



/// \brief fast access site data structure -- formatted for use in the
/// lhood calculation loop
///
struct site_data_fastlup : public site_data_fastlup_core, private uncopyable{

  typedef site_data_fastlup_core base_t;

  site_data_fastlup()
    : base_t(),n_assigned_data_sets(0),total_count(0),
      n_groups(0),group_label(),data(0),_is_fixed(false) {}

  ~site_data_fastlup() { clear(); };

  void init(unsigned init_len,
            unsigned init_n_orgs,
            unsigned init_n_assigned_data_sets);

  /// sort and setup for fast indexing
  void fix();

  /// \brief ads mask is true for sites excluded in each ads. (I guess
  /// it's actually an inverse mask)
  ///
  /// this is used to make assigned data set lhood calculations more efficient
  ///
  const std::vector<bool>&
  ads_site_mask(const unsigned a) const {
    return _ads_site_mask[a];
  }

private:
  void clear(bool is_clear_base=true);

public:
  // think about listing cat_counts in single full list, with sites,
  // so that this info can be looped over in one pass with no zeros
  unsigned n_assigned_data_sets;
  unsigned total_count;
  unsigned n_groups;
  std::vector<std::string> group_label;
  site_data_fastlup_cell* data;

private:
  std::vector<std::vector<bool> > _ads_site_mask;
  bool _is_fixed;
};


#endif

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

// $Id: site_data.h 1193 2008-03-29 03:14:54Z ctsa $

/// \file
///
/// cat and group info is stored in sites
///


#ifndef __SITE_DATA_H
#define __SITE_DATA_H

#include "audit_info.h"
#include "site_data_code.h"
#include "name_id_lup.h"
#include "site_code.h"
#include "util/bio/bioseq_util.h"

#include <iosfwd>
#include <map>
#include <string>
#include <utility>


typedef std::pair<const site_data_code,unsigned> site_count_pair;

#ifdef NO_BOOST_POOL_ALLOCATOR
#ifdef __GNUC__
#include <ext/pool_allocator.h>
typedef __gnu_cxx::__pool_alloc<site_count_pair> site_count_alloc;
#else
typedef std::allocator<site_count_pair> site_count_alloc;
#endif
#else
#include <boost/pool/pool_alloc.hpp>
typedef boost::fast_pool_allocator<site_count_pair> site_count_alloc;
#endif

typedef std::map<site_data_code,unsigned,std::less<site_data_code>, site_count_alloc> site_count_data_type;

typedef std::map<site_code,site_count_data_type> site_count_map;




/// \brief holds processed sequence site data in an exchangeable form
///
struct site_data {

  site_data() : sm(SITE_MODEL::NONE), _ai() {}

  explicit
  site_data(const audit_info& ai) : sm(SITE_MODEL::NONE), _ai(ai) {}

  unsigned n_taxa() const { return taxid.size(); }

  unsigned data_class_size() const { return data_class_labels.size(); }

  // quick debugging
  void peek(std::ostream& os) const;

  //
  void store_state(std::ostream& os) const;

  //
  void load_state(std::istream& is,
                  const bool is_read_group_info = true);

  void clear() {
    data_class_labels.clear();
    taxid.clear();
    site_count.clear();
    sm=SITE_MODEL::NONE;
    group_label.clear();
  }

  //////////////////////////////////
  name_id_lup data_class_labels;
  name_id_lup taxid;
  site_count_map site_count;
  SITE_MODEL::index_t sm;
  std::vector<std::string> group_label;

private:
  const audit_info _ai;
};


#endif

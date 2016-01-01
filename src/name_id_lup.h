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

// $Id: name_id_lup.h 945 2007-10-23 23:44:18Z ctsa $

/// \file

#ifndef __NAME_ID_LUP_H
#define __NAME_ID_LUP_H

#include <map>
#include <string>
#include <vector>


/// \brief maps between id nums and string labels
struct name_id_lup{

  name_id_lup() : idmap(), strmap() {}

  /// \brief add label if not present, and return label id
  unsigned assignid(const std::string& s);

  /// \brief get id of pre-existing label
  unsigned getid(const std::string& s) const;

  /// \brief does label exist?
  bool testid(const std::string& s) const;

  /// \brief get pre-existing label
  const std::string& getstr(const unsigned id) const;

  unsigned size() const { return strmap.size(); }

  void clear() { idmap.clear(); strmap.clear(); }

private:
  std::map<std::string,unsigned> idmap;
  std::vector<std::string> strmap;
};


#endif

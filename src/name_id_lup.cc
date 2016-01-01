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

// $Id: name_id_lup.cc 945 2007-10-23 23:44:18Z ctsa $

/// \file

#include "name_id_lup.h"
#include "substk_exception.h"

#include <map>
#include <string>
#include <vector>


unsigned
name_id_lup::
assignid(const std::string& s) {
  const std::map<std::string,unsigned>::const_iterator a(idmap.find(s));
  if(a==idmap.end()){
    const unsigned id(strmap.size());
    idmap[s]=id;
    strmap.push_back(s);
    return id;
  } else {
    return a->second;
  }
}



bool
name_id_lup::
testid(const std::string& s) const {
   const std::map<std::string,unsigned>::const_iterator a(idmap.find(s));
   return !(a==idmap.end());
}



unsigned
name_id_lup::
getid(const std::string& s) const {
  const std::map<std::string,unsigned>::const_iterator a(idmap.find(s));
  if(a==idmap.end()) throw substk_exception("name_id_lup.getid(): can't find name");
  return a->second;
}



const std::string&
name_id_lup::
getstr(const unsigned id) const {
  if(id>=strmap.size()) throw substk_exception("name_id_lup.getid(): invalid id");
  return strmap[id];
}

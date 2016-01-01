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

// $Id: string_util.h 942 2007-10-23 22:39:35Z ctsa $

/// \file

#include <string>

/// strip leading and trailing w/s
///
inline
void strip(std::string& s) {
  static const char ws[] = " \r\t\n";

  if(s.empty()) return;
  std::string::size_type start(s.find_first_not_of(ws));
  std::string::size_type end(s.find_last_not_of(ws));
  s = s.substr(start,(end-start+1));
}


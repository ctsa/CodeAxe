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

// $Id: audit_info.h 755 2007-08-17 02:31:46Z ctsa $

/// \file

#ifndef __AUDIT_INFO_H
#define __AUDIT_INFO_H

#include <string>
#include <iosfwd>


/// \brief holds optional auditing info about the program run
///
struct audit_info {

  audit_info() : seed(-1) {}

  void store_state(std::ostream& os) const;

  void skip_state(std::istream& is) const;

  std::string cmdline;
  std::string bininfo;
  long int seed;
};


#endif

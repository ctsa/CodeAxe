// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// SubsTK : phylogenetic analysis and simulation library
//
//   http://www.phrap.org
//
//
// Copyright (C) 2007 Christopher T Saunders (ctsa@u.washington.edu)
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

// $Id: io_util.h 1108 2008-01-25 03:04:10Z ctsa $

/// \file

#ifndef __IO_UTIL_H
#define __IO_UTIL_H

#include <iosfwd>

/// \brief check non-empty input fstream
void
check_nonempty_istream(std::istream& is,
                       const char* const stream_label = "UNKNOWN");

/// \brief clear out ws, stop after newline
void
clear_whitespace_line(std::istream& is);

/// \brief clear to newline
void
clear_line(std::istream& is);

/// \brief read int from is (cleaner code look for temporaries)
int
read_int(std::istream& is);


struct iliner {

  iliner(const char* const s) : section_id(s) {}

  void
  advance(std::istream& is,
          const char* const second = 0) const;

private:
  const char* const section_id;
};

#endif

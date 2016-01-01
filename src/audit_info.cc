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

// $Id: audit_info.cc 1134 2008-02-11 22:18:55Z ctsa $

/// \file

#include "audit_info.h"
#include "util/general/io_util.h"

#include <ctime>

#include <iostream>


const char section_id[] = "audit_info";
const char end_label[] = "END";
const iliner il(section_id);


void
audit_info::
store_state(std::ostream& os) const {
  os << section_id << " cmdline " << cmdline << "\n";
  os << section_id << " binary " << bininfo << "\n";
  os << section_id << " random_seed " << seed << "\n";

  const time_t result(time(0));
  os << section_id << " write_time " << asctime(localtime(&result));
  os << section_id << " " << end_label << "\n";
}



static
void
skip_advance(std::istream& is,
             const char* const second){
  il.advance(is,second);
  clear_line(is);
}



void
audit_info::
skip_state(std::istream& is) const {
  skip_advance(is,"cmdline");
  skip_advance(is,"binary");
  skip_advance(is,"random_seed");
  skip_advance(is,"write_time");
  skip_advance(is,end_label);
}

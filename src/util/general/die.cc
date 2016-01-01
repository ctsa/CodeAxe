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

// $Id: die.cc 892 2007-10-02 18:45:30Z ctsa $

/// \file

#include "die.h"
#include "log.h"

#include <cstdlib>

#include <ostream>


void
die(const char * msg){
  log_os << "ERROR:: " << msg << std::endl;
  abort();
}


void
pass_away(const char * msg){
  log_os << "ERROR:: " << msg << std::endl;
  exit(EXIT_FAILURE);
}


void
warning(const char * msg){
  log_os << "WARNING:: " << msg << std::endl;
}


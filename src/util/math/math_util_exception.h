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

// $Id: math_util_exception.h 931 2007-10-22 19:07:34Z ctsa $

/// \file

#ifndef __MATH_UTIL_EXCEPTION_H
#define __MATH_UTIL_EXCEPTION_H

#include "../general/die.h"

#include <exception>
#include <string>


class math_util_exception : public std::exception {
 public:
  math_util_exception(const char* s) : message(s) {
#ifdef KILL_EXCEPTIONS
    die(s);
#endif
  }
  ~math_util_exception() throw() {};

  const char* what() const throw() { return message.c_str(); }

  std::string message;
};


#endif


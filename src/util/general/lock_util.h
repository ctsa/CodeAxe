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

// $Id: lock_util.h 1026 2007-11-26 18:54:47Z ctsa $

/// \file

#ifndef __LOCK_UTIL_H
#define __LOCK_UTIL_H

#include "uncopyable.h"


/// \brief returns file descriptor
///
int
set_file_lock(const char* file);


///
void
unset_file_lock(const int fd);

/// \brief simple file locking: lifetime of class == lifetime of lock
///
struct file_locker : private uncopyable {
  file_locker(const char* filename)
    : fd(set_file_lock(filename)) {}

  ~file_locker() { unset_file_lock(fd); }

private:
  int fd;
};

#endif

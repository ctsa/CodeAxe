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

// $Id: lock_util.cc 731 2007-08-12 17:38:19Z ctsa $

/// \file

#include "die.h"
#include "lock_util.h"

#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
//#include <sys/file.h>

int
set_file_lock(const char* file){

#ifdef DARWIN_HACK
  struct flock fl = { 0, 0, getpid(), F_WRLCK, SEEK_SET };
#else
  struct flock fl = { F_WRLCK, SEEK_SET, 0, 0, getpid() };
#endif

  int fd;
  if((fd = open(file, O_WRONLY|O_CREAT,S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH)) == -1){
    die("file open failed");
  }

  if(fcntl(fd, F_SETLKW, &fl) == -1) { perror("fcntl: "); die("file lock failed"); }  /* F_GETLK, F_SETLK, F_SETLKW */
  //if(lockf(fd,F_LOCK,0) == -1) { die("file lock failed"); }

  return fd;
}


void
unset_file_lock(const int fd){

#ifdef DARWIN_HACK
  struct flock fl = { 0, 0, getpid(), F_UNLCK, SEEK_SET };
#else
  struct flock fl = { F_UNLCK, SEEK_SET, 0, 0, getpid() };
#endif

  if(fcntl(fd, F_SETLK, &fl) == -1) { perror("fcntl: "); die("file unlock failed"); }
    //if(lockf(fd,F_ULOCK,0) == -1) { die("file unlock failed"); }

  if(close(fd) == -1) { die("file close failed"); }
}

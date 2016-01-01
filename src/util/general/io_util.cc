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

// $Id: io_util.cc 1108 2008-01-25 03:04:10Z ctsa $

/// \file

#include "io_util.h"
#include "log.h"
#include "metatags.h"

#include <cstdlib>

#include <iostream>
#include <string>

void
check_nonempty_istream(std::istream& is,
                       const char* stream_label){

  // rudimentary stream check
  if( ! is ) {
    log_os << "ERROR:: invalid input stream: " << stream_label << "\n";
    exit(EXIT_FAILURE);
  }
  is.get();
  if( is.eof() ){
    log_os << "ERROR:: input stream is empty: " << stream_label << "\n";
    exit(EXIT_FAILURE);
  }
  is.unget();
}


// clear out ws, stop after newline
void
clear_whitespace_line(std::istream& is){
  while(is){
    const int c(is.get());
    if(c==' ' || c== '\t') continue;

    if(c!='\n') is.unget();
    break;
  }
}

// clear to newline
void
clear_line(std::istream& is){
  while(is){
    const int c(is.get());

    if(c=='\n') { break;}
  }
}


int
read_int(std::istream& is) {
  int i;
  is >> i;
  return i;
}



static
void
iliner_die(std::istream& is,
           const char* const sid,
           const char* const tmpstr)  NORETURN_TAG;

static
void
iliner_die(std::istream& is,
           const char* const sid,
           const char* const tmpstr){
  
  log_os << "ERROR:: -- Invalid model file format. --\n"
         << "ERROR:: expected_section: " << sid << "\n"
         << "ERROR:: section: " << tmpstr << "\n"
         << "\n";
  log_os << is.rdbuf();
  abort();
}



void
iliner::
advance(std::istream& is,
        const char* const second) const {

  std::string tmpstr;
  is >> tmpstr;
  if(tmpstr != section_id) iliner_die(is,section_id,tmpstr.c_str());
  is >> tmpstr;
  if(second && tmpstr != second) iliner_die(is,second,tmpstr.c_str());
}

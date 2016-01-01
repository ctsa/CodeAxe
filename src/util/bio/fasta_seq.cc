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

// $Id: fasta_seq.cc 585 2007-06-12 20:16:09Z ctsa $

/// \file 

#include "fasta_seq.h"

#include <iterator>
#include <iostream>

using namespace std;


istream& operator>>(istream& is,FastaSeq& f){
  getline(is,f.header);

  f.seq.clear();
  istream_iterator<char> i(is),i_end;

  while(!(i == i_end) && *i != '>'){
    f.seq.push_back(*i);
    ++i;
  }
  if(is.eof() && f.seq.size()) { is.clear(); }
  is.unget();
  return is;
}


ostream& operator<<(ostream& os, const FastaSeq& f){
  static const unsigned line_length(70);

  os << f.header;
  vector<char>::const_iterator i,i_begin=f.seq.begin(),i_end=f.seq.end();
  for(i=i_begin;i!=i_end;++i){
    if(! ((i-i_begin)%line_length)) os << "\n";
    os << *i;
  }
  os << "\n";
  return os;
}

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

// $Id: fasta_seq.h 585 2007-06-12 20:16:09Z ctsa $

/// \file 

#ifndef __FASTA_SEQ_H
#define __FASTA_SEQ_H

#include <iosfwd>
#include <string>
#include <vector>

/// \brief simple fasta object with slow but easy to use i/o
///
struct FastaSeq {
  FastaSeq() : header(), seq() {}
  FastaSeq(const std::string& header_, const std::vector<char>& seq_)
    : header(header_), seq(seq_) {}
  std::string header;
  std::vector<char> seq;
};

std::istream& operator>>(std::istream& is,FastaSeq& f);

std::ostream& operator<<(std::ostream& os,const FastaSeq& f);

#endif

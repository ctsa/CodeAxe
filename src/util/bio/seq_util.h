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

// $Id: seq_util.h 734 2007-08-13 23:53:46Z ctsa $

/// \file 
///
/// \brief some char based bioseq utils, enumeration versions in
/// bioseq_util are preferred for speed
///

#ifndef __SEQ_UTIL_H
#define __SEQ_UTIL_H

#include <algorithm>
#include <functional>
#include <iomanip>
#include <ostream>
#include <iterator>
#include <map>
#include <string>


template <typename SeqIter>
void count_seq_members(SeqIter first,
                       SeqIter last,
                       std::map<typename std::iterator_traits<SeqIter>::value_type,int>& counts){

  for(;first!=last;++first) counts[*first]++;
}

template <typename SeqIter>
void count_seq_triplets(SeqIter first,
                        SeqIter last,
                        std::map<std::string,int>& c){
  typedef typename std::iterator_traits<SeqIter>::value_type InputT;
//  boost::function_requires<boost::ConvertibleConcept<InputT,char> >();

  InputT c0,c1,c2;
  c0 = c1 = c2 = 0;
  for(unsigned count=0;first!=last;++first,++count){
    switch(count%3){
     case 0:
      c0 = *first; break;
     case 1:
      c1 = *first; break;
     case 2:
      c2 = *first;
      char a[] = { c0,c1,c2,0 };
      c[std::string(a)]++;
    }
  }
}


struct base_complement : public std::unary_function<char,char> {
  char operator()(char c){
    switch(c){
     case 'A': return 'T';
     case 'C': return 'G';
     case 'G': return 'C';
     case 'T': return 'A';
     default : return 'N';
    }
  }
};


template <typename SeqIter>
void reverse_complement(SeqIter first,
                        SeqIter last){
  typedef typename std::iterator_traits<SeqIter>::value_type InputT;
//  boost::function_requires<boost::ConvertibleConcept<InputT,char> >();

  std::transform(first,last,first,base_complement());
  std::reverse(first,last);
}



template <typename InfoType>
void write_codon_info(std::ostream& os, std::map<std::string,InfoType>& t, unsigned w=6){

  // make table
  const char print_order[5] = "TCAG";
  std::string codon = "XXX";
  for(int row=0;row<16;++row){
    codon[0] = print_order[row/4];
    codon[2] = print_order[row%4];
    for(int col=0;col<4;++col){
      codon[1] = print_order[col];
      os << " " << std::setw(w) << t[codon];
    }
    os << "\n";
    if( row%4 == 3 ) os << "\n";
  }
}


// true if c points to a 4fold codon
bool is_4fold_codon(const char c0,const char c1);
bool is_4fold_codon(const char* c);


// ambiguous character codon translator
char codon_trans(const char* c);

// literal type1 codon translator
char codon_trans_known(const char* c);

#endif

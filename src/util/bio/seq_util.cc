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

// $Id: seq_util.cc 585 2007-06-12 20:16:09Z ctsa $

/// \file 

#include "seq_util.h"

#include <cstring>

#include <map>
#include <utility>


using namespace std;


const unsigned CODON_BASE_SIZE(3);


struct ltstr {
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};


struct ltstr_codon {
  bool operator()(const char* s1, const char* s2) const
  {
    return strncmp(s1, s2,CODON_BASE_SIZE) < 0;
  }
};


// return true for 4fold codons
//
bool
is_4fold_codon(const char c0,const char c1){
  if ( c0 == 'C' || c0 == 'G' ) {
    if ( c1 == 'T' || c1 == 'C' || c1 == 'G' ) return true;
  } else if ( c0 == 'T' || c0 == 'A' ) {
    if ( c1 == 'C' ) return true;
  }
  return false;
}

bool is_4fold_codon(const char* c) {
  return is_4fold_codon(c[0],c[1]);
}




// general ambiguous codon translator
//
char codon_trans(const char*c){

  typedef map<char,const char*> amg_map;

  const static pair<char,const char*> tab_tmp[] = {
    make_pair('A',"A"),
    make_pair('C',"C"),
    make_pair('G',"G"),
    make_pair('T',"T"),
    make_pair('K',"GT"),
    make_pair('M',"AC"),
    make_pair('R',"AG"),
    make_pair('S',"GC"),
    make_pair('W',"AT"),
    make_pair('Y',"CT"),
    make_pair('B',"CGT"),
    make_pair('D',"AGT"),
    make_pair('H',"ACT"),
    make_pair('V',"ACG"),
    make_pair('N',"ACGT"),
    make_pair('X',"ACGT")};

  const static int tsize(sizeof(tab_tmp)/sizeof(tab_tmp[0]));

  const static amg_map amg(tab_tmp,tab_tmp+tsize);

  if ( strncmp(c,"---",CODON_BASE_SIZE) == 0 ) return '-';

  char result(codon_trans_known(c));

  if(result != 'X') return result;

  amg_map::const_iterator s0 = amg.find(c[0]);
  amg_map::const_iterator s1 = amg.find(c[1]);
  amg_map::const_iterator s2 = amg.find(c[2]);

  if(s0 == amg.end() || s1 == amg.end() || s2 == amg.end()) return 'X';


  char old_result(0);
  char codon[CODON_BASE_SIZE];
  for(const char* c0=s0->second;*c0;++c0){
    codon[0] = *c0;
    for(const char* c1=s1->second;*c1;++c1){
      codon[1] = *c1;
      for(const char* c2=s2->second;*c2;++c2){
        codon[2] = *c2;
        result = codon_trans_known(codon);
        if ( old_result && result != old_result ) return 'X';
        old_result = result;
      }
    }
  }

  return result;
}




// literal type1 codon translator
//
char codon_trans_known(const char*c){
  typedef map<const char*,char,ltstr_codon> lup_map;

  const static pair<const char*,char> tab_tmp[] = {
    make_pair("TTT",'F'),make_pair("TCT",'S'),make_pair("TAT",'Y'),make_pair("TGT",'C'),
    make_pair("TTC",'F'),make_pair("TCC",'S'),make_pair("TAC",'Y'),make_pair("TGC",'C'),
    make_pair("TTA",'L'),make_pair("TCA",'S'),make_pair("TAA",'*'),make_pair("TGA",'*'),
    make_pair("TTG",'L'),make_pair("TCG",'S'),make_pair("TAG",'*'),make_pair("TGG",'W'),

    make_pair("CTT",'L'),make_pair("CCT",'P'),make_pair("CAT",'H'),make_pair("CGT",'R'),
    make_pair("CTC",'L'),make_pair("CCC",'P'),make_pair("CAC",'H'),make_pair("CGC",'R'),
    make_pair("CTA",'L'),make_pair("CCA",'P'),make_pair("CAA",'Q'),make_pair("CGA",'R'),
    make_pair("CTG",'L'),make_pair("CCG",'P'),make_pair("CAG",'Q'),make_pair("CGG",'R'),

    make_pair("ATT",'I'),make_pair("ACT",'T'),make_pair("AAT",'N'),make_pair("AGT",'S'),
    make_pair("ATC",'I'),make_pair("ACC",'T'),make_pair("AAC",'N'),make_pair("AGC",'S'),
    make_pair("ATA",'I'),make_pair("ACA",'T'),make_pair("AAA",'K'),make_pair("AGA",'R'),
    make_pair("ATG",'M'),make_pair("ACG",'T'),make_pair("AAG",'K'),make_pair("AGG",'R'),

    make_pair("GTT",'V'),make_pair("GCT",'A'),make_pair("GAT",'D'),make_pair("GGT",'G'),
    make_pair("GTC",'V'),make_pair("GCC",'A'),make_pair("GAC",'D'),make_pair("GGC",'G'),
    make_pair("GTA",'V'),make_pair("GCA",'A'),make_pair("GAA",'E'),make_pair("GGA",'G'),
    make_pair("GTG",'V'),make_pair("GCG",'A'),make_pair("GAG",'E'),make_pair("GGG",'G')};

  const static int tsize(sizeof(tab_tmp)/sizeof(tab_tmp[0]));

  const static lup_map codon_table1(tab_tmp,tab_tmp+tsize);

  if ( strncmp(c,"---",CODON_BASE_SIZE) == 0 ) return '-';

  // gcc4 has some type of mystery inlining error right here:
  const lup_map::const_iterator s(codon_table1.find(c));
  if( s != codon_table1.end() ){
    return s->second;
  } else {
    return 'X';
  }
}

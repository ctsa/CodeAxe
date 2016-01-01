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

// $Id: bioseq_util.cc 675 2007-07-10 20:47:03Z ctsa $

/// \file 

#include "bioseq_util.h"

#include <cassert>

#include <algorithm>
#include <ostream>



const char NUC::syms[] = {'A','C','G','T','N'};

const char AA::syms[] = {'A','C','D','E','F','G','H','I','K','L',
                         'M','N','P','Q','R','S','T','V','W','Y','*','X'};

const char NSAA::syms[] = {'A','C','D','E','F','G','H','I','K','L',
                           'M','N','P','Q','R','S','T','V','W','Y','X'};




namespace DINUC{
  std::ostream&
  operator<<(std::ostream& os,
             const print& p){
    NUC::index_t nx = DINUC::decode_nx(p._c);
    NUC::index_t n0 = DINUC::decode_n0(p._c);

    os << NUC::syms[nx] << "+" << NUC::syms[n0];
    return os;
  }
}


namespace CODON{
  std::ostream&
  operator<<(std::ostream& os,
             const print& p){
    NUC::index_t n[BASE_SIZE];
    CODON::decode(n,p._c);
    os << NUC::syms[n[0]] << NUC::syms[n[1]] << NUC::syms[n[2]];
    if(p._is_with_aa) os << "(" << AA::syms[codon_trans_known(p._c)] << ")";
    return os;
  }
}


namespace NSCODON{
  std::ostream&
  operator<<(std::ostream& os,
             const print& p){
    NUC::index_t n[CODON::BASE_SIZE];
    NSCODON::decode(n,p._c);
    os << NUC::syms[n[0]] << NUC::syms[n[1]] << NUC::syms[n[2]];
    if(p._is_with_aa) os << "(" << NSAA::syms[codon_trans_known(p._c)] << ")";
    return os;
  }
}


namespace SITE_MODEL{
  std::ostream&
  operator<<(std::ostream& os,
             const print& p){

    switch(p._sm){
    case AA:
      os << AA::syms[p._state];
      break;
    case NSAA:
      os << NSAA::syms[p._state];
      break;
    default:
      NUC::index_t n[MAX_BASE_SIZE];
      unsigned cpos[MAX_BASE_SIZE];

      const unsigned n_bases(base_size(p._sm));
      decode_nuc(p._sm,p._state,n);
      codon_position(p._sm,cpos);
      for(unsigned i(0);i<n_bases;++i){
        if(i && cpos[i] < cpos[i-1]){ os << "+"; }
        os << NUC::syms[n[i]];
      }
    }
    if(is_coding(p._sm) && p._is_tcode){
      os << "(" << AA::syms[codon_trans_known(decode_codon(p._sm,p._state))] << ")";
    }
    return os;
  }
}



bool
is_4fold_codon(const NUC::index_t n0,\
               const NUC::index_t n1){
  if ( n0 == NUC::C || n0 == NUC::G ) {
    if ( n1 == NUC::T || n1 == NUC::C || n1 == NUC::G ) return true;
  } else if ( n0 == NUC::T || n0 == NUC::A ) {
    if ( n1 == NUC::C ) return true;
  }
  return false;
}

bool
is_4fold_codon(const NUC::index_t* n) {
  return is_4fold_codon(n[0],n[1]);
}


#ifdef DARWIN_HACK
// gcc in os x 10.4 can't seem to call the constructor before allowing
// calls to codon_table->translate. that's just scary.
//
static bool ctinit(false);
#endif


struct codon_table {

  codon_table(){ init(ct); }

  AA::index_t translate(const CODON::index_t c) {
#ifdef DARWIN_HACK
    if(!ctinit) { init(ct); ctinit=true; }
#endif
    return ct[c];
  }

private:
  static
  void init(AA::index_t ct[]);

  AA::index_t ct[CODON::SIZE];
};



void
codon_table::
init(AA::index_t ct[]){

  using namespace CODON;
  using namespace AA;

#ifdef FLAT_CODON_TABLE
  for(unsigned i(0);i<CODON::SIZE;++i){ ct[i] = A; }
#else
  ct[TTT]=F; ct[TCT]=S; ct[TAT]=Y;    ct[TGT]=C;
  ct[TTC]=F; ct[TCC]=S; ct[TAC]=Y;    ct[TGC]=C;
  ct[TTA]=L; ct[TCA]=S; ct[TAA]=STOP; ct[TGA]=STOP;
  ct[TTG]=L; ct[TCG]=S; ct[TAG]=STOP; ct[TGG]=W;
  
  ct[CTT]=L; ct[CCT]=P; ct[CAT]=H; ct[CGT]=R;
  ct[CTC]=L; ct[CCC]=P; ct[CAC]=H; ct[CGC]=R;
  ct[CTA]=L; ct[CCA]=P; ct[CAA]=Q; ct[CGA]=R;
  ct[CTG]=L; ct[CCG]=P; ct[CAG]=Q; ct[CGG]=R;
  
  ct[ATT]=I; ct[ACT]=T; ct[AAT]=N; ct[AGT]=S;
  ct[ATC]=I; ct[ACC]=T; ct[AAC]=N; ct[AGC]=S;
  ct[ATA]=I; ct[ACA]=T; ct[AAA]=K; ct[AGA]=R;
  ct[ATG]=M; ct[ACG]=T; ct[AAG]=K; ct[AGG]=R;
  
  ct[GTT]=V; ct[GCT]=A; ct[GAT]=D; ct[GGT]=G;
  ct[GTC]=V; ct[GCC]=A; ct[GAC]=D; ct[GGC]=G;
  ct[GTA]=V; ct[GCA]=A; ct[GAA]=E; ct[GGA]=G;
  ct[GTG]=V; ct[GCG]=A; ct[GAG]=E; ct[GGG]=G;
#endif
}



codon_table ctab;



/// high speed type1 codon table -- does not tolerate substitutions
///
AA::index_t
codon_trans_known(const CODON::index_t c){
  return ctab.translate(c);
}



void
aa_reverse_trans(CODON::index_t c[CODON::MAX_CODON_PER_AA],
                 unsigned& n_codons,
                 const AA::index_t a){
  n_codons = 0;
  for(unsigned i(0);i<CODON::SIZE;++i){
    const AA::index_t aa(ctab.translate(static_cast<CODON::index_t>(i)));
    if(a == aa){
      c[n_codons++] = static_cast<CODON::index_t>(i);
    }
  }
}


NSAA::index_t
codon_trans_known(NSCODON::index_t c){
  return NSAA::convert_aa_to_nsaa(codon_trans_known(NSCODON::convert_nscodon_to_codon(c)));
}


void
aa_reverse_trans(NSCODON::index_t c[CODON::MAX_CODON_PER_AA],
                 unsigned& n_codons,
                 const NSAA::index_t a){
  CODON::index_t cc[CODON::MAX_CODON_PER_AA];
  aa_reverse_trans(cc,n_codons,convert_nsaa_to_aa(a));
  for(unsigned i(0);i<n_codons;++i){
    c[i] = NSCODON::convert_codon_to_nscodon(cc[i]);
  }
}


/// given two aa's, enumerate all possible single nuc change codon paths
/// between them::
void
get_codon_paths(std::vector<nscodon_pair>& vc,
                const NSAA::index_t aa1,
                const NSAA::index_t aa2){
  vc.clear();

  if( aa1 == aa2 ) return;

  unsigned nc1,nc2;
  NSCODON::index_t cl1[CODON::MAX_CODON_PER_AA];
  NSCODON::index_t cl2[CODON::MAX_CODON_PER_AA];

  aa_reverse_trans(cl1,nc1,aa1);
  aa_reverse_trans(cl2,nc2,aa2);

  // determine which codon pairs are separated by a single nuc:
  for(unsigned cl1i(0);cl1i<nc1;++cl1i){
    NSCODON::index_t c1 = cl1[cl1i];
    NUC::index_t c1n1,c1n2,c1n3;
    NSCODON::decode(c1n1,c1n2,c1n3,static_cast<NSCODON::index_t>(c1));

    for(unsigned cl2i(0);cl2i<nc2;++cl2i){
      NSCODON::index_t c2 = cl2[cl2i];
      assert(c1!=c2);

      NUC::index_t c2n1,c2n2,c2n3;
      NSCODON::decode(c2n1,c2n2,c2n3,static_cast<NSCODON::index_t>(c2));

      unsigned nsubs(0);
      NUC::index_t c1ndiff(NUC::N);
      NUC::index_t c2ndiff(NUC::N);

      if( c1n1 != c2n1 ) {nsubs++; c1ndiff=c1n1; c2ndiff=c2n1;}
      if( c1n2 != c2n2 ) {nsubs++; c1ndiff=c1n2; c2ndiff=c2n2;}
      if( c1n3 != c2n3 ) {nsubs++; c1ndiff=c1n3; c2ndiff=c2n3;}

      if(nsubs == 1){
        vc.push_back(nscodon_pair(c1,c2));
      }
    }
  }
}


void
get_first_nuc_diff(NUC::index_t n1,
                   NUC::index_t n2,
                   unsigned n,
                   const NSCODON::index_t c1,
                   const NSCODON::index_t c2){
  n1 = NUC::N;
  n2 = NUC::N;
  n = 0;

  if(c1==c2) return;

  NUC::index_t c1n[3];
  NSCODON::decode(c1n[0],c1n[1],c1n[2],static_cast<NSCODON::index_t>(c1));

  NUC::index_t c2n[3];
  NSCODON::decode(c2n[0],c2n[1],c2n[2],static_cast<NSCODON::index_t>(c2));

  for(unsigned i(0);i<3;++i){
    if(c1n[i]!=c2n[i]){
      n=i;
      n1=c1n[i];
      n2=c2n[i];
      return;
    }
  }
}


// matrix of which aa's reach other aa's by single mutations
void
make_reachable_aa_matrix(bool is_reachable[NSAA::SIZE*NSAA::SIZE]){

  for(unsigned i(0);i<NSAA::SIZE*NSAA::SIZE;++i)
    is_reachable[i] = false;

  // build list of excluded transitions, same for reversible/non-reversible models
  for(unsigned c1(0);c1<NSCODON::SIZE;++c1){
    NUC::index_t c1n1,c1n2,c1n3;
    NSCODON::decode(c1n1,c1n2,c1n3,static_cast<NSCODON::index_t>(c1));
    const NSAA::index_t c1aa =
      static_cast<NSAA::index_t>(codon_trans_known(NSCODON::convert_nscodon_to_codon(static_cast<NSCODON::index_t>(c1))));

    for(unsigned c2(c1+1);c2<NSCODON::SIZE;++c2){
      if( c1==c2 ) continue;
      NUC::index_t c2n1,c2n2,c2n3;
      NSCODON::decode(c2n1,c2n2,c2n3,static_cast<NSCODON::index_t>(c2));
      const NSAA::index_t c2aa =
        static_cast<NSAA::index_t>(codon_trans_known(NSCODON::convert_nscodon_to_codon(static_cast<NSCODON::index_t>(c2))));
      if( c1aa == c2aa ) continue;

      unsigned nsubs(0);

      if( c1n1 != c2n1 ) nsubs++;
      if( c1n2 != c2n2 ) nsubs++;
      if( c1n3 != c2n3 ) nsubs++;

      if( nsubs == 1 ) {
        is_reachable[c2aa+c1aa*NSAA::SIZE] = true;
        is_reachable[c1aa+c2aa*NSAA::SIZE] = true;
      }
    }
  }
}

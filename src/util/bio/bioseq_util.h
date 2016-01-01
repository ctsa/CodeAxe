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

// $Id: bioseq_util.h 786 2007-09-06 00:45:00Z ctsa $

/// \file 
/// \brief definitions of biological state types 
///

#ifndef __BIOSEQ_UTIL_H
#define __BIOSEQ_UTIL_H

#include "../general/die.h"

#include <functional>
#include <iosfwd>
#include <vector>


namespace NUC {
  enum index_t { A,C,G,T,N };
  static const index_t SIZE(N);
  static const index_t AMBIG(N);

  extern const char syms[SIZE+1];

  inline
  bool
  is_transition(const index_t n1,const index_t n2){
    return ((n1==A && n2==G) || (n1==G && n2==A) ||
            (n1==C && n2==T) || (n1==T && n2==C));
  }

  inline
  index_t
  remap(const char a){
    for(unsigned i=0;i<SIZE;++i) {
      if( a == syms[i] ) return static_cast<index_t>(i);
    }
    return N;
  }
}

// (nx,n0) nuc pair
namespace DINUC {
  typedef unsigned index_t;
  enum size_t { SIZE=NUC::SIZE*NUC::SIZE };
  static const index_t AMBIG(SIZE);

  inline
  NUC::index_t
  decode_nx(index_t d) { return static_cast<NUC::index_t>(d/NUC::SIZE); }

  inline
  NUC::index_t
  decode_n0(index_t d) { return static_cast<NUC::index_t>(d%NUC::SIZE); }

  inline
  index_t
  encode(NUC::index_t nx,
         NUC::index_t n0) { return n0+nx*NUC::SIZE; }


  // stream dinuc printer::
  class print {
    index_t _c;
  public:
    print(const index_t c = 0)
      : _c(c) {}

    friend
    std::ostream&
    operator<<(std::ostream& os,
               const print& p);
  };
}

namespace AA {
  enum index_t { A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,STOP,X };
  enum size_t { SIZE=X };
  extern const char syms[SIZE+1];
  static const size_t AMBIG(SIZE);
}

/// non-stop aa to go with non-stop codon
namespace NSAA {
  enum index_t { A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,X };
  enum size_t { SIZE=X };
  extern const char syms[SIZE+1];
  static const size_t AMBIG(SIZE);

  inline
  index_t
  convert_aa_to_nsaa(AA::index_t a){
    if(a == AA::STOP || a == AA::X) return X;
    else return static_cast<index_t>(a);
  }

  inline
  AA::index_t
  convert_nsaa_to_aa(index_t a){
    if(a == X) return AA::X;
    else return static_cast<AA::index_t>(a);
  }
}

namespace CODON {
  static const unsigned BASE_SIZE(3);
  static const unsigned MAX_CODON_PER_AA(6);

  enum index_t {
    AAA,CAA,GAA,TAA, ACA,CCA,GCA,TCA, AGA,CGA,GGA,TGA, ATA,CTA,GTA,TTA,
    AAC,CAC,GAC,TAC, ACC,CCC,GCC,TCC, AGC,CGC,GGC,TGC, ATC,CTC,GTC,TTC,
    AAG,CAG,GAG,TAG, ACG,CCG,GCG,TCG, AGG,CGG,GGG,TGG, ATG,CTG,GTG,TTG,
    AAT,CAT,GAT,TAT, ACT,CCT,GCT,TCT, AGT,CGT,GGT,TGT, ATT,CTT,GTT,TTT, NNN };

  static const index_t SIZE(NNN);
  static const index_t AMBIG(NNN);


  inline
  index_t
  encode(const NUC::index_t n1,
         const NUC::index_t n2,
         const NUC::index_t n3){

    if(n1==NUC::N || n2==NUC::N || n3==NUC::N) {
      return NNN;
    } else {
      return static_cast<index_t>((n3*NUC::SIZE+n2)*NUC::SIZE+n1);
    }
  }

  template <typename RandomAccessIterator>
  index_t
  encode(RandomAccessIterator n){
    return encode(n[0],n[1],n[2]);
  }

  inline
  void
  decode(NUC::index_t& n1,
         NUC::index_t& n2,
         NUC::index_t& n3,
         const index_t c){
    if(c==NNN) {
      n1 = n2 = n3 = NUC::N;
      return;
    }
    n1 = static_cast<NUC::index_t>(c%NUC::SIZE);

    const unsigned tmp(c/NUC::SIZE);
    n2 = static_cast<NUC::index_t>(tmp%NUC::SIZE);
    n3 = static_cast<NUC::index_t>(tmp/NUC::SIZE);
  }

  inline
  void
  decode(NUC::index_t* n,
         const index_t c){
    decode(n[0],n[1],n[2],c);
  }

#if 0
  inline
  NUC::index_t
  decode_nuc3(const index_t c){
    return static_cast<NUC::index_t>(c/(NUC::SIZE*NUC::SIZE));
  }
#endif

  // stream codon printer::
  class print {
    index_t _c;
    bool _is_with_aa;
  public:
    print(const index_t c = NNN,
          const bool is_with_aa = false)
      : _c(c), _is_with_aa(is_with_aa) {}

    friend
    std::ostream&
    operator<<(std::ostream& os,
               const print& p);
  };
}


// make a new name for CODON namespace to refer to noncoding 3-mers:
namespace TRINUC { using namespace CODON; }


/// non-stop codon space
///
namespace NSCODON {
#ifdef FLAT_CODON_TABLE
  enum index_t {
    AAA,CAA,GAA,TAA, ACA,CCA,GCA,TCA, AGA,CGA,GGA,TGA, ATA,CTA,GTA,TTA,
    AAC,CAC,GAC,TAC, ACC,CCC,GCC,TCC, AGC,CGC,GGC,TGC, ATC,CTC,GTC,TTC,
    AAG,CAG,GAG,TAG, ACG,CCG,GCG,TCG, AGG,CGG,GGG,TGG, ATG,CTG,GTG,TTG,
    AAT,CAT,GAT,TAT, ACT,CCT,GCT,TCT, AGT,CGT,GGT,TGT, ATT,CTT,GTT,TTT, NNN };
#else
  enum index_t {
    AAA,CAA,GAA,     ACA,CCA,GCA,TCA, AGA,CGA,GGA,     ATA,CTA,GTA,TTA,
    AAC,CAC,GAC,TAC, ACC,CCC,GCC,TCC, AGC,CGC,GGC,TGC, ATC,CTC,GTC,TTC,
    AAG,CAG,GAG,     ACG,CCG,GCG,TCG, AGG,CGG,GGG,TGG, ATG,CTG,GTG,TTG,
    AAT,CAT,GAT,TAT, ACT,CCT,GCT,TCT, AGT,CGT,GGT,TGT, ATT,CTT,GTT,TTT, NNN };
#endif
  static const index_t SIZE(NNN);
  static const index_t AMBIG(NNN);

#ifdef FLAT_CODON_TABLE
  inline
  index_t
  convert_codon_to_nscodon(const CODON::index_t c){

    if(c == CODON::NNN) return NNN;

    unsigned precodon(static_cast<unsigned>(c));
    return static_cast<index_t>(precodon);
  }

  inline
  CODON::index_t
  convert_nscodon_to_codon(const index_t c){

    if(c==NNN) return CODON::NNN;

    unsigned precodon(static_cast<unsigned>(c));
    return static_cast<CODON::index_t>(precodon);
  }
#else
  inline
  index_t
  convert_codon_to_nscodon(const CODON::index_t c){

    if(c == CODON::NNN) return NNN;

    unsigned precodon(static_cast<unsigned>(c));

    if      (precodon> CODON::TAG) { precodon -= 3; }
    else if (precodon==CODON::TAG) { return NNN; }
    else if (precodon> CODON::TGA) { precodon -= 2; }
    else if (precodon==CODON::TGA) { return NNN; }
    else if (precodon> CODON::TAA) { precodon -= 1; }
    else if (precodon==CODON::TAA) { return NNN; }
    return static_cast<index_t>(precodon);
  }

  inline
  CODON::index_t
  convert_nscodon_to_codon(const index_t c){

    if(c==NNN) return CODON::NNN;

    unsigned precodon(static_cast<unsigned>(c));

    if      (precodon+3> CODON::TAG) { precodon += 3;}
    else if (precodon+2> CODON::TGA) { precodon += 2;}
    else if (precodon+1> CODON::TAA) { precodon += 1;}

    return static_cast<CODON::index_t>(precodon);
  }
#endif

  inline
  index_t
  encode(const NUC::index_t n1,
         const NUC::index_t n2,
         const NUC::index_t n3){
    return convert_codon_to_nscodon(CODON::encode(n1,n2,n3));
  }

  template <typename RandomAccessIterator>
  index_t
  encode(RandomAccessIterator n){
    return encode(n[0],n[1],n[2]);
  }

  inline
  void
  decode(NUC::index_t& n1,
         NUC::index_t& n2,
         NUC::index_t& n3,
         index_t c){
    CODON::decode(n1,n2,n3,convert_nscodon_to_codon(c));
  }

  inline
  void
  decode(NUC::index_t* n,
         const index_t c){
    decode(n[0],n[1],n[2],c);
  }

#if 0
  inline
  NUC::index_t
  decode_nuc3(const index_t c){
    return CODON::decode_nuc3(convert_nscodon_to_codon(c));
  }
#endif

  // stream codon printer::
  class print {
    index_t _c;
    bool _is_with_aa;
  public:
    print(const index_t c = NNN,
          const bool is_with_aa = false)
      : _c(c), _is_with_aa(is_with_aa) {}

    friend
    std::ostream&
    operator<<(std::ostream& os,
               const print& p);
  };
}




// nsc4 is nscodon plus either the preceding or following nucleotide
// (nx)
namespace NSC4 {
  typedef unsigned index_t;

  static const index_t SIZE(NSCODON::SIZE*NUC::SIZE);
  static const index_t AMBIG(SIZE);

  inline
  NUC::index_t
  decode_nuc(const index_t i){
    if(i==AMBIG){
      return NUC::N;
    } else {
      return static_cast<NUC::index_t>(i/NSCODON::SIZE);
    }
  }

  inline
  NSCODON::index_t
  decode_nscodon(const index_t i){
    if(i==AMBIG){
      return NSCODON::NNN;
    } else {
      return static_cast<NSCODON::index_t>(i%NSCODON::SIZE);
    }
  }

  inline
  void
  decode(NUC::index_t& n,
         NSCODON::index_t& c,
         const index_t i){

    n = decode_nuc(i);
    c = decode_nscodon(i);
  }

  inline
  index_t
  encode(const NUC::index_t n,
         const NSCODON::index_t c){
    if(n==NUC::N || c==NSCODON::NNN){
      return AMBIG;
    }else{
      return n*NSCODON::SIZE+c;
    }
  }

  /// for a given n, provide a continuous id offset for c:
  inline
  unsigned
  get_nscodon_offset(const NUC::index_t n){
    return encode(n,static_cast<NSCODON::index_t>(0));
  }
}




//
namespace NSC5 {
  typedef unsigned index_t;

  static const index_t SIZE(NSC4::SIZE*NUC::SIZE);
  static const index_t AMBIG(SIZE);

  inline
  NUC::index_t
  decode_nuc1(const index_t i){
    if(i==AMBIG){
      return NUC::N;
    } else {
      return static_cast<NUC::index_t>((i/NSCODON::SIZE)%NUC::SIZE);
    }
  }

  inline
  NUC::index_t
  decode_nuc2(const index_t i){
    if(i==AMBIG){
      return NUC::N;
    } else {
      return static_cast<NUC::index_t>((i/NSCODON::SIZE)/NUC::SIZE);
    }
  }

  inline
  NSCODON::index_t
  decode_nscodon(const index_t i){
    if(i==AMBIG){
      return NSCODON::NNN;
    } else {
      return static_cast<NSCODON::index_t>(i%NSCODON::SIZE);
    }
  }

  inline
  void
  decode(NUC::index_t& n1,
         NUC::index_t& n2,
         NSCODON::index_t& c,
         const index_t i){
    n1 = decode_nuc1(i);
    n2 = decode_nuc2(i);
    c = decode_nscodon(i);
  }

  inline
  index_t
  encode(const NUC::index_t n1,
         const NUC::index_t n2,
         const NSCODON::index_t c){
    if(n1==NUC::N || n2==NUC::N || c==NSCODON::NNN){
      return AMBIG;
    } else {
      return ((n2*NUC::SIZE+n1)*NSCODON::SIZE+c);
    }
  }

  /// for a given (n_pre,post), provide a continous id offset for c:
  inline
  unsigned
  get_nscodon_offset(const NUC::index_t n1,
                     const NUC::index_t n2){
    return encode(n1,n2,static_cast<NSCODON::index_t>(0));
  }
}



/// programatic access to the bioseq site types:
namespace SITE_MODEL {
  enum index_t { NONE, NUC, DINUC, TRINUC, CODON, NSCODON, AA, NSAA, NSC4PRE, NSC4POST, NSC5, NSC5PRE, NSC5POST, BINARY };

  inline
  unsigned
  base_size(const index_t i){
    switch(i){
    case NUC:      return 1;
    case DINUC:    return 2;
    case TRINUC:   return 3;
    case CODON:
    case NSCODON:  return CODON::BASE_SIZE;
    case NSC4PRE:
    case NSC4POST: return CODON::BASE_SIZE+1;
    case NSC5PRE:
    case NSC5POST:
    case NSC5:     return CODON::BASE_SIZE+2;
    default:       return 0;
    }
  }

  inline
  unsigned
  state_size(const index_t i){
    switch(i){
    case BINARY:   return 2;
    case NUC:      return NUC::SIZE;
    case DINUC:    return NUC::SIZE*NUC::SIZE;
    case TRINUC:   return CODON::SIZE;
    case CODON:    return CODON::SIZE;
    case NSCODON:  return NSCODON::SIZE;
    case NSC4PRE:
    case NSC4POST: return NSC4::SIZE;
    case NSC5PRE:
    case NSC5POST:
    case NSC5:     return NSC5::SIZE;
    default:       return 0;
    }
  }

  inline
  unsigned
  is_coding(const index_t i){
    switch(i){
    case CODON:
    case NSCODON:
    case NSC4PRE:
    case NSC4POST:
    case NSC5PRE:
    case NSC5POST:
    case NSC5:     return true;
    default:       return false;
    }
  }

  inline
  unsigned
  ambig_state(const index_t i){
    switch(i){
    case NUC:      return NUC::AMBIG;
    case DINUC:    return DINUC::AMBIG;
    case TRINUC:   return CODON::AMBIG;
    case CODON:    return CODON::AMBIG;
    case NSCODON:  return NSCODON::AMBIG;
    case NSC4PRE:
    case NSC4POST: return NSC4::AMBIG;
    case NSC5PRE:
    case NSC5POST:
    case NSC5:     return NSC5::AMBIG;
    case AA:       return AA::AMBIG;
    case NSAA:     return NSAA::AMBIG;
    default:       die("unknown ambig state for this site model");
    }
  }

  // relieve client code from heap-allocing base_size arrays...
  enum { MAX_BASE_SIZE=CODON::BASE_SIZE+2 };

  inline
  void
  codon_position(const index_t i,
                 unsigned* n){

    const unsigned b(base_size(i));

    unsigned start_pos(0);
    switch(i){
    case CODON:
    case NSCODON:
    case NSC4POST:
    case NSC5POST: start_pos=0; break;
    case NSC4PRE:
    case NSC5:     start_pos=2; break;
    case NSC5PRE:  start_pos=1; break;
    default:
      for(unsigned j(0);j<b;++j){ n[j]=0; }
      return;
    }

    for(unsigned j(0);j<b;++j){ n[j]=(j+start_pos)%CODON::BASE_SIZE; }
  }

  inline
  void
  is_partial_codon(const index_t i,
                   bool* n){

    const unsigned b(base_size(i));
    for(unsigned j(0);j<b;++j){ n[j]=false; }

    switch(i){
    case CODON:
    case NSCODON: break;
    case NSC4PRE:  n[0]=true; break;
    case NSC4POST: n[3]=true; break;
    case NSC5:     n[0]=true; n[4]=true; break;
    case NSC5PRE:  n[0]=true; n[1]=true; break;
    case NSC5POST: n[3]=true; n[4]=true; break;
    default: die("invalid index to is_partial_codon");
    }
  }


  template <typename RandomAccessIterator>
  unsigned
  encode_nuc(const index_t i,
             RandomAccessIterator nuc){

    NSCODON::index_t c;

    switch(i){
    case NUC:    return nuc[0];
    case DINUC:  return DINUC::encode(nuc[0],nuc[1]);
    case TRINUC: return TRINUC::encode(nuc);
    case CODON:  return CODON::encode(nuc);
    case NSCODON: return NSCODON::encode(nuc);
    case NSC4PRE:
      c = NSCODON::encode(nuc+1);
      return NSC4::encode(nuc[0],c);
    case NSC4POST:
      c = NSCODON::encode(nuc);
      return NSC4::encode(nuc[3],c);
    case NSC5PRE:
      c = NSCODON::encode(nuc+2);
      return NSC5::encode(nuc[0],nuc[1],c);
    case NSC5:
      c = NSCODON::encode(nuc+1);
      return NSC5::encode(nuc[0],nuc[4],c);
    case NSC5POST:
      c = NSCODON::encode(nuc);
      return NSC5::encode(nuc[3],nuc[4],c);
    default: die("invalid index to decode_nuc");
    }
  }

  inline
  void
  decode_nuc(const index_t i,
             const unsigned state,
             NUC::index_t* n){

    NSCODON::index_t c;

    switch(i){
    case NUC:
      n[0] = static_cast<NUC::index_t>(state);
      break;
    case DINUC:
      n[0] = DINUC::decode_nx(static_cast<DINUC::index_t>(state));
      n[1] = DINUC::decode_n0(static_cast<DINUC::index_t>(state));
      break;
    case TRINUC:
      TRINUC::decode(n,static_cast<TRINUC::index_t>(state));
      break;
    case CODON:
      CODON::decode(n,static_cast<CODON::index_t>(state));
      break;
    case NSCODON:
      NSCODON::decode(n,static_cast<NSCODON::index_t>(state));
      break;
    case NSC4PRE:
      NSC4::decode(n[0],c,static_cast<NSC4::index_t>(state));
      NSCODON::decode(n+1,c);
      break;
    case NSC4POST:
      NSC4::decode(n[3],c,static_cast<NSC4::index_t>(state));
      NSCODON::decode(n,c);
      break;
    case NSC5PRE:
      NSC5::decode(n[0],n[1],c,static_cast<NSC5::index_t>(state));
      NSCODON::decode(n+2,c);
      break;
    case NSC5:
      NSC5::decode(n[0],n[4],c,static_cast<NSC5::index_t>(state));
      NSCODON::decode(n+1,c);
      break;
    case NSC5POST:
      NSC5::decode(n[3],n[4],c,static_cast<NSC5::index_t>(state));
      NSCODON::decode(n,c);
      break;
    default: die("invalid index to decode_nuc");
    }
  }

  inline
  NSCODON::index_t
  decode_nscodon(const index_t i,
                 const unsigned state){
    switch(i){
    case CODON:    return NSCODON::convert_codon_to_nscodon(static_cast<CODON::index_t>(state));
    case NSCODON:  return static_cast<NSCODON::index_t>(state);
    case NSC4PRE:
    case NSC4POST: return NSC4::decode_nscodon(static_cast<NSC4::index_t>(state));
    case NSC5PRE:
    case NSC5POST:
    case NSC5:     return NSC5::decode_nscodon(static_cast<NSC5::index_t>(state));
    default: die("invalid index to decode_nscodon");
    }
  }

  inline
  CODON::index_t
  decode_codon(const index_t i,
               const unsigned state){
    switch(i){
    case CODON:    return static_cast<CODON::index_t>(state);
    case NSCODON:
    case NSC4PRE:
    case NSC4POST:
    case NSC5PRE:
    case NSC5POST:
    case NSC5:     return NSCODON::convert_nscodon_to_codon(decode_nscodon(i,state));
    default:  die("invalid index to decode_codon");
    }
  }

  inline
  const char*
  label(const index_t i){

    switch(i){
    case NUC:      return "nuc";
    case DINUC:    return "dinuc";
    case TRINUC:   return "trinuc";
    case CODON:    return "codon";
    case NSCODON:  return "nscodon";
    case NSC4PRE:  return "nsc4-pre";
    case NSC4POST: return "nsc4-post";
    case NSC5:     return "nsc5-center";
    case NSC5PRE:  return "nsc5-pre";
    case NSC5POST: return "nsc5-post";
    default: die("invalid index to decode_nscodon");
    }
  }

  // stream site printer
  struct print {
    print(const index_t sm = NONE,
          const unsigned state = 0,
          const bool is_translate_coding = false)
      : _sm(sm), _state(state), _is_tcode(is_translate_coding) {}

    friend
    std::ostream&
    operator<<(std::ostream& os,
               const print& p);

  private:
    index_t _sm;
    unsigned _state;
    bool _is_tcode;
  };
}



// how a module would be defined in doxygen, but don't want it here anymore:
//
//  /// \defgroup NucUtils Nucleotide Utilities
//  ///
//  //@{

struct nuc_index_complement : public std::unary_function<NUC::index_t,NUC::index_t> {
  NUC::index_t operator()(NUC::index_t n){
    using namespace NUC;
    switch(n){
     case A:  return T;
     case C:  return G;
     case G:  return C;
     case T:  return A;
     default: return N;
    }
  }
};

//  //@}


bool
is_4fold_codon(const NUC::index_t n0,
               const NUC::index_t n1);

bool
is_4fold_codon(const NUC::index_t* n);


AA::index_t
codon_trans_known(CODON::index_t c);

void
aa_reverse_trans(CODON::index_t i[CODON::MAX_CODON_PER_AA],
                 unsigned& n_codons,
                 const AA::index_t a);

NSAA::index_t
codon_trans_known(NSCODON::index_t c);


void
aa_reverse_trans(NSCODON::index_t i[CODON::MAX_CODON_PER_AA],
                 unsigned& n_codons,
                 const NSAA::index_t a);


// matrix of which aa's reach other aa's by single mutations
//
void
make_reachable_aa_matrix(bool is_reachable[NSAA::SIZE*NSAA::SIZE]);


struct nscodon_pair {
  nscodon_pair(NSCODON::index_t _c1,
               NSCODON::index_t _c2) : c1(_c1), c2(_c2) {}
  NSCODON::index_t c1;
  NSCODON::index_t c2;
};


/// given two aa's, enumerate all possible single nuc change codon paths
/// between them::
void
get_codon_paths(std::vector<nscodon_pair>& vc,
                const NSAA::index_t aa1,
                const NSAA::index_t aa2);


/// given two codons, list the first differing nucleotide and position,
/// return N's and position == 0 otherwise.
///
void
get_first_nuc_diff(NUC::index_t n1,
                   NUC::index_t n2,
                   unsigned n,
                   const NSCODON::index_t c1,
                   const NSCODON::index_t c2);

#endif

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

// $Id: bioseq_util_pdf_conversions.hh 1056 2007-12-09 22:32:55Z ctsa $

/// \file 

#include "../math/prob_util.h"
#include "../math/prob_util_io.h"
#include "../general/log.h"



template <typename RAI,
          typename RAI2>
void
bio_pdistro_convert(const BIO_PDISTRO::index_t from_index,
                    const BIO_PDISTRO::index_t to_index,
                    const RAI f,
                    RAI2 t){

  if(from_index==to_index){
    std::copy(f,f+BIO_PDISTRO::param_size(from_index),t);
    return;
  }

  switch(from_index){
  case BIO_PDISTRO::NUC:
    switch(to_index){
    case BIO_PDISTRO::DINUC:    nuc_pdf_2_dinuc_pdf(f,t);   return;
    case BIO_PDISTRO::TRINUC:   nuc_pdf_2_trinuc_pdf(f,t);  return;
    case BIO_PDISTRO::CODON:    nuc_pdf_2_codon_pdf(f,t);   return;
    case BIO_PDISTRO::NSCODON:  nuc_pdf_2_nscodon_pdf(f,t); return;
    case BIO_PDISTRO::NSC4PRE:  nuc_pdf_2_nsc4_pdf(f,t);    return;
    case BIO_PDISTRO::NSC4POST: nuc_pdf_2_nsc4_pdf(f,t);    return;
    default: break;
    }
    break;
  case BIO_PDISTRO::NUC_CODON_POS:
    switch(to_index){
    case BIO_PDISTRO::NSCODON:  nuc_pdf_all_pos_2_nscodon_pdf(f,t);   return;
    case BIO_PDISTRO::NSC4PRE:  nuc_pdf_all_pos_2_nsc4_pre_pdf(f,t);  return;
    case BIO_PDISTRO::NSC4POST: nuc_pdf_all_pos_2_nsc4_post_pdf(f,t); return;
    default: break;
    }
    break;
  case BIO_PDISTRO::DINUC:
    switch(to_index){
    case BIO_PDISTRO::NUC: dinuc_pdf_2_nuc_pdf_pos1(f,t); return;
    default: break;
    }
    break;
  case BIO_PDISTRO::TRINUC:
    switch(to_index){
    case BIO_PDISTRO::NUC: trinuc_pdf_2_nuc_pdf_pos(f,2,t); return;
    default: break;
    }
    break;
  case BIO_PDISTRO::NSCODON:
    switch(to_index){
    case BIO_PDISTRO::NUC:           nscodon_pdf_2_nuc_pdf(f,t);         return;
    case BIO_PDISTRO::NUC_CODON_POS: nscodon_pdf_2_nuc_pdf_all_pos(f,t); return;
    case BIO_PDISTRO::NSAA:          nscodon_pdf_2_nsaa_pdf(f,t);        return;
    case BIO_PDISTRO::NSC4PRE:       nscodon_pdf_2_nsc4_pre_pdf(f,t);    return;
    case BIO_PDISTRO::NSC4POST:      nscodon_pdf_2_nsc4_post_pdf(f,t);   return;
    case BIO_PDISTRO::NSC5:          nscodon_pdf_2_nsc5_pdf(f,t);        return;
    default: break;
    }
    break;
  case BIO_PDISTRO::NSC4PRE:
  case BIO_PDISTRO::NSC4POST:
    switch(to_index){
    case BIO_PDISTRO::NUC:           nsc4_pdf_2_nuc_pdf(f,t);         return;
    case BIO_PDISTRO::NUC_CODON_POS: nsc4_pdf_2_nuc_pdf_all_pos(f,t); return;
    case BIO_PDISTRO::NSCODON:       nsc4_pdf_2_nscodon_pdf(f,t);     return;
    default: break;
    }
    break;
  case BIO_PDISTRO::NSC5:
  case BIO_PDISTRO::NSC5PRE:
  case BIO_PDISTRO::NSC5POST:
    switch(to_index){
    case BIO_PDISTRO::NSCODON: nsc5_pdf_2_nscodon_pdf(f,t); return;
    default: break;
    }
    break;
  default: break;
  }

  die("invalid distro conversion combo requested in bio_pdistro_convert");
}




template <typename RAI>
void
#ifdef DEBUG
bioseq_util_pdistro_check(const RAI pdf,
                          const unsigned size){
  static const double BIOSEQ_PTOL(1.e-5);
  pdistro_check(pdf,size,BIOSEQ_PTOL);
}
#else
bioseq_util_pdistro_check(const RAI,
                          const unsigned){}
#endif




template <typename RAI,
          typename RAI2>
void
nuc_pdf_2_dinuc_pdf(const RAI nuc_pdf,
                    RAI2 dinuc_pdf){

  for(int i(0);i<DINUC::SIZE;++i){
    int nx = DINUC::decode_nx(i);
    int n0 = DINUC::decode_n0(i);
    *(dinuc_pdf+i) = *(nuc_pdf+nx) * *(nuc_pdf+n0);
  }

  bioseq_util_pdistro_check(dinuc_pdf,DINUC::SIZE);
}




template <typename RAI,
          typename RAI2>
void
nuc_pdf_2_codon_pdf(const RAI nuc_pdf,
                    RAI2 codon_pdf){

  NUC::index_t nuc[CODON::BASE_SIZE];
  for(int i(0);i<CODON::SIZE;++i){
    CODON::decode(nuc,static_cast<CODON::index_t>(i));
    *(codon_pdf+i) = *(nuc_pdf+nuc[0]) * *(nuc_pdf+nuc[1]) * *(nuc_pdf+nuc[2]);
  }

  bioseq_util_pdistro_check(codon_pdf,CODON::SIZE);
}




template <typename RAI,
          typename RAI2>
void
nuc_pdf_2_nscodon_pdf(const RAI nuc_pdf, //[NUC::SIZE],
                      RAI2 nscodon_pdf){ //[NSCODON::SIZE]){

  NUC::index_t n1,n2,n3;

  for(unsigned c(0);c<NSCODON::SIZE;++c){
    NSCODON::decode(n1,n2,n3,static_cast<NSCODON::index_t>(c));
    *(nscodon_pdf+c) = *(nuc_pdf+n1) * *(nuc_pdf+n2) * *(nuc_pdf+n3);
  }
  pdistro_norm(nscodon_pdf,nscodon_pdf+NSCODON::SIZE);

  bioseq_util_pdistro_check(nscodon_pdf,NSCODON::SIZE);
}




template <typename RAI,
          typename RAI2>
void
nuc_pdf_2_nsc4_pdf(const RAI nuc_pdf,
                   RAI2 nsc4_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  T nscodon_pdf[NSCODON::SIZE];
  nuc_pdf_2_nscodon_pdf(nuc_pdf,nscodon_pdf);
  nscodon_pdf_nuc_pdf_2_nsc4_pdf(nscodon_pdf,nuc_pdf,nsc4_pdf);
}




template <typename RAI,
          typename RAI2>
void
nuc_pdf_all_pos_2_codon_pdf(const RAI nuc_pdf,
                            RAI2 codon_pdf){

  NUC::index_t n1,n2,n3;
  for(unsigned c(0);c<CODON::SIZE;++c){
    CODON::decode(n1,n2,n3,static_cast<CODON::index_t>(c));
    *(codon_pdf+c) = *(nuc_pdf+n1) * *(nuc_pdf+n2+NUC::SIZE) * *(nuc_pdf+n3+NUC::SIZE*2);
  }
  pdistro_norm(codon_pdf,codon_pdf+CODON::SIZE);

  bioseq_util_pdistro_check(codon_pdf,CODON::SIZE);
}




template <typename RAI,
          typename RAI2>
void
nuc_pdf_all_pos_2_nscodon_pdf(const RAI nuc_pdf,
                              RAI2 nscodon_pdf){

  // because of stop codons, the nscodon_pdf has nuc pos priors,
  // this is so that we can recover something closer to the generative all_nuc distro
  // from nscodon_pdf later, although this is still not exact:
  //
  typedef typename std::iterator_traits<RAI>::value_type T;

  T nuc_prior_pdf[NUC::SIZE*CODON::BASE_SIZE];
  pdistro_unif(nscodon_pdf,nscodon_pdf+NSCODON::SIZE);
  nscodon_pdf_2_nuc_pdf_all_pos(nscodon_pdf,nuc_prior_pdf);
  for(unsigned i(0);i<CODON::BASE_SIZE;++i){
    for(unsigned j(0);j<NUC::SIZE;++j){
      nuc_prior_pdf[NUC::SIZE*i+j] = *(nuc_pdf+NUC::SIZE*i+j) / *(nuc_prior_pdf+NUC::SIZE*i+j);
    }
    pdistro_norm(nuc_prior_pdf+NUC::SIZE*i,nuc_prior_pdf+NUC::SIZE*(i+1));
  }

  NUC::index_t n1,n2,n3;
  for(unsigned c(0);c<NSCODON::SIZE;++c){
    NSCODON::decode(n1,n2,n3,static_cast<NSCODON::index_t>(c));
    //nscodon_pdf[c] = nuc_pdf[n1]*nuc_pdf[n2+NUC::SIZE]*nuc_pdf[n3+NUC::SIZE*2];
    *(nscodon_pdf+c) = nuc_prior_pdf[n1]*nuc_prior_pdf[n2+NUC::SIZE]*nuc_prior_pdf[n3+NUC::SIZE*2];
  }
  pdistro_norm(nscodon_pdf,nscodon_pdf+NSCODON::SIZE);

  bioseq_util_pdistro_check(nscodon_pdf,NSCODON::SIZE);
}




template <typename RAI,typename RAI2>
void
nuc_pdf_all_pos_2_nsc4_pdf_pos(const RAI nuc_pdf,
                               const unsigned pos,
                               RAI2 nsc4_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  T nscodon_pdf[NSCODON::SIZE];
  nuc_pdf_all_pos_2_nscodon_pdf(nuc_pdf,nscodon_pdf);
  nscodon_pdf_nuc_pdf_2_nsc4_pdf(nscodon_pdf,nuc_pdf+pos*NUC::SIZE,nsc4_pdf);
}




template <typename RAI,typename RAI2>
void
nuc_pdf_all_pos_2_nsc4_pre_pdf(const RAI nuc_pdf,
                               RAI2 nsc4_pdf){
  nuc_pdf_all_pos_2_nsc4_pdf_pos(nuc_pdf,2,nsc4_pdf);
}




template <typename RAI,typename RAI2>
void
nuc_pdf_all_pos_2_nsc4_post_pdf(const RAI nuc_pdf,
                                RAI2 nsc4_pdf){
  nuc_pdf_all_pos_2_nsc4_pdf_pos(nuc_pdf,0,nsc4_pdf);
}




template <typename RAI,typename RAI2>
void
dinuc_pdf_2_nuc_pdf(const RAI dinuc_pdf,
                    RAI2 nuc_pdf){

  for(unsigned i(0);i<NUC::SIZE;++i) *(nuc_pdf+i) = 0.;

  for(int i(0);i<DINUC::SIZE;++i){
    const int nx(DINUC::decode_nx(i));
    const int n0(DINUC::decode_n0(i));
    *(nuc_pdf+nx) += *(dinuc_pdf+i)/2.;
    *(nuc_pdf+n0) += *(dinuc_pdf+i)/2.;
  }

  bioseq_util_pdistro_check(nuc_pdf,NUC::SIZE);
}




template <typename RAI,typename RAI2>
void
dinuc_pdf_2_nuc_pdf_pos1(const RAI dinuc_pdf,
                         RAI2 nuc_pdf){

  for(unsigned i(0);i<NUC::SIZE;++i) *(nuc_pdf+i) = 0.;

  for(int i(0);i<DINUC::SIZE;++i){
    const int n0(DINUC::decode_n0(i));
    *(nuc_pdf+n0) += *(dinuc_pdf+i);
  }

  bioseq_util_pdistro_check(nuc_pdf,NUC::SIZE);
}




template <typename RAI,typename FloatType>
void
dinuc_pdf_2_conditioned_dinuc_pdf(const RAI dinuc_pdf,
                                  const unsigned condition_position,
                                  FloatType conditioned_dinuc_pdf[NUC::SIZE][NUC::SIZE]){

  const unsigned p(condition_position);

  NUC::index_t nuc[2];
  for(unsigned n_cond(0);n_cond<NUC::SIZE;++n_cond){
    nuc[(0+p)%2] = static_cast<NUC::index_t>(n_cond);
    FloatType sum(0.);
    for(unsigned n(0);n<NUC::SIZE;++n){
      nuc[(1+p)%2] = static_cast<NUC::index_t>(n);
      sum += dinuc_pdf[DINUC::encode(nuc[0],nuc[1])];
    }

    for(unsigned n(0);n<NUC::SIZE;++n){
      nuc[(1+p)%2] = static_cast<NUC::index_t>(n);
      conditioned_dinuc_pdf[n_cond][n] = dinuc_pdf[DINUC::encode(nuc[0],nuc[1])]/sum;
    }
  }
}




template <typename RAI>
void
dinuc_pdf_2_conditioned_dinuc_pdf_inplace(const unsigned condition_position,
                                          RAI dinuc_pdf){

  typedef typename std::iterator_traits<RAI>::value_type T;

  const unsigned p(condition_position);

  NUC::index_t nuc[2];
  for(unsigned n_cond(0);n_cond<NUC::SIZE;++n_cond){
    nuc[(0+p)%2] = static_cast<NUC::index_t>(n_cond);
    T sum(0.);
    for(unsigned n(0);n<NUC::SIZE;++n){
      nuc[(1+p)%2] = static_cast<NUC::index_t>(n);
      sum += *(dinuc_pdf+DINUC::encode(nuc[0],nuc[1]));
    }

    for(unsigned n(0);n<NUC::SIZE;++n){
      nuc[(1+p)%2] = static_cast<NUC::index_t>(n);
      const DINUC::index_t di(DINUC::encode(nuc[0],nuc[1]));
      *(dinuc_pdf+di) /= sum;
    }
  }
}




template <typename RAI,typename RAI2>
void
dinuc_pdf_2_conditioned_dinuc_pdf(const RAI dinuc_pdf,
                                  const unsigned condition_position,
                                  RAI2 conditioned_dinuc_pdf){

  std::copy(dinuc_pdf,dinuc_pdf+DINUC::SIZE,conditioned_dinuc_pdf);
  dinuc_pdf_2_conditioned_dinuc_pdf_inplace(condition_position,conditioned_dinuc_pdf);
}




template <typename RAI,typename RAI2>
void
codon_pdf_2_nuc_pdf_pos(const RAI codon_pdf,
                        const unsigned pos,
                        RAI2 nuc_pdf){

  if(pos >= CODON::BASE_SIZE) {
    log_os << "Invalid codon position: " << pos << "\n";
    abort();
  }

  for(unsigned i(0);i<NUC::SIZE;++i) *(nuc_pdf+i) = 0.;

  NUC::index_t n[CODON::BASE_SIZE];
  for(unsigned i(0);i<CODON::SIZE;++i){
    CODON::decode(n,static_cast<CODON::index_t>(i));
    *(nuc_pdf+n[pos]) += *(codon_pdf+i);
  }

  bioseq_util_pdistro_check(nuc_pdf,NUC::SIZE);
}




template <typename RAI,typename RAI2>
void
codon_pdf_2_nuc_pdf_all_pos(const RAI codon_pdf,
                            RAI2 nuc_pdf){

  for(unsigned i(0);i<NUC::SIZE*CODON::BASE_SIZE;++i) *(nuc_pdf+i) = 0.;

  NUC::index_t n[CODON::BASE_SIZE];
  for(unsigned i(0);i<CODON::SIZE;++i){
    CODON::decode(n,static_cast<CODON::index_t>(i));
    for(unsigned pos(0);pos<CODON::BASE_SIZE;++pos){
      *(nuc_pdf+NUC::SIZE*pos+n[pos]) += *(codon_pdf+i);
    }
  }

  for(unsigned i(0);i<CODON::BASE_SIZE;++i){
    bioseq_util_pdistro_check(nuc_pdf+NUC::SIZE*i,NUC::SIZE);
  }
}




template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nuc_pdf(const RAI nscodon_pdf,
                      RAI2 nuc_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  for(unsigned i(0);i<NUC::SIZE;++i) *(nuc_pdf+i) = 0.;

  NUC::index_t nuc[CODON::BASE_SIZE];
  for(unsigned i(0);i<NSCODON::SIZE;++i){
    const T& pi(*(nscodon_pdf+i));
    NSCODON::decode(nuc,static_cast<NSCODON::index_t>(i));
    for(unsigned j(0);j<CODON::BASE_SIZE;++j){ *(nuc_pdf+nuc[j]) += pi; }
  }

  pdistro_norm(nuc_pdf,nuc_pdf+NUC::SIZE);
  bioseq_util_pdistro_check(nuc_pdf,NUC::SIZE);
}




template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nuc_pdf_pos(const RAI nscodon_pdf,
                          const unsigned pos,
                          RAI2 nuc_pdf){

  if(pos >= CODON::BASE_SIZE) {
    log_os << "Invalid codon position: " << pos << "\n";
    abort();
  }

  for(unsigned i(0);i<NUC::SIZE;++i) *(nuc_pdf+i) = 0.;

  NUC::index_t n[CODON::BASE_SIZE];
  for(unsigned i(0);i<NSCODON::SIZE;++i){
    NSCODON::decode(n,static_cast<NSCODON::index_t>(i));
    *(nuc_pdf+n[pos]) += *(nscodon_pdf+i);
  }

  bioseq_util_pdistro_check(nuc_pdf,NUC::SIZE);
}



template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nuc_pdf_all_pos(const RAI nscodon_pdf,
                              RAI2 nuc_pdf){

  for(unsigned i(0);i<NUC::SIZE*CODON::BASE_SIZE;++i) *(nuc_pdf+i) = 0.;

  NUC::index_t n[CODON::BASE_SIZE];
  for(unsigned i(0);i<NSCODON::SIZE;++i){
    NSCODON::decode(n,static_cast<NSCODON::index_t>(i));
    for(unsigned pos(0);pos<CODON::BASE_SIZE;++pos){
      *(nuc_pdf+NUC::SIZE*pos+n[pos]) += *(nscodon_pdf+i);
    }
  }

  for(unsigned i(0);i<CODON::BASE_SIZE;++i){
    bioseq_util_pdistro_check(nuc_pdf+NUC::SIZE*i,NUC::SIZE);
  }
}



template <typename FloatType>
void
nscodon_pdf_2_dinuc_pdfs_conditioned_on_pos(const FloatType nscodon_pdf[NSCODON::SIZE],
                                            const int pos,
                                            FloatType conditioned_nscodon_pdf[NUC::SIZE][DINUC::SIZE]){

  int pos1((pos+1)%CODON::BASE_SIZE);
  int pos2((pos+2)%CODON::BASE_SIZE);

  if(pos1>pos2) std::swap(pos1,pos2);

  NUC::index_t nuc[CODON::BASE_SIZE];
  for(int n(0);n<NUC::SIZE;++n){
    nuc[pos] = static_cast<NUC::index_t>(n);
    for(int d(0);d<DINUC::SIZE;++d) { conditioned_nscodon_pdf[n][d] = 0.; }
    FloatType sum(0.);
    for(int d(0);d<DINUC::SIZE;++d){
      nuc[pos1] = DINUC::decode_nx(d);
      nuc[pos2] = DINUC::decode_n0(d);
      const NSCODON::index_t c(NSCODON::encode(nuc));
      if(c==NSCODON::NNN) continue;
      sum += nscodon_pdf[c];
    }
    for(int d(0);d<DINUC::SIZE;++d){
      nuc[pos1] = DINUC::decode_nx(d);
      nuc[pos2] = DINUC::decode_n0(d);
      const NSCODON::index_t c(NSCODON::encode(nuc));
      if(c==NSCODON::NNN) {
        conditioned_nscodon_pdf[n][d] = 0.;
      } else {
        conditioned_nscodon_pdf[n][d] = nscodon_pdf[c]/sum;
      }
    }
  }

  for(int n(0);n<NUC::SIZE;++n){
    bioseq_util_pdistro_check(conditioned_nscodon_pdf[n],DINUC::SIZE);
  }
}




template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nsaa_pdf(const RAI nscodon_pdf,
                       RAI2 nsaa_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  for(unsigned i(0);i<NSAA::SIZE;++i) *(nsaa_pdf+i) = 0.;

  for(unsigned i(0);i<NSCODON::SIZE;++i){
    const T& pi(*(nscodon_pdf+i));

    NSAA::index_t aa;
    aa = codon_trans_known(static_cast<NSCODON::index_t>(i));
    *(nsaa_pdf+aa) += pi;
  }

  bioseq_util_pdistro_check(nsaa_pdf,NSAA::SIZE);
}




template <typename RAI,typename RAI2,typename RAI3>
void
nscodon_pdf_nuc_pdf_2_nsc4_pdf(const RAI nscodon_pdf,
                               const RAI2 nuc_pdf,
                               RAI3 nsc4_pdf){

  for(unsigned i(0);i<NSC4::SIZE;++i){
    const NSC4::index_t fi(static_cast<NSC4::index_t>(i));
    const NUC::index_t ni(NSC4::decode_nuc(fi));
    const NSCODON::index_t ci(NSC4::decode_nscodon(fi));
    *(nsc4_pdf+i) = *(nscodon_pdf+ci) * *(nuc_pdf+ni);
  }

  bioseq_util_pdistro_check(nsc4_pdf,NSC4::SIZE);
}




template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nsc4_pdf_allnuc(const RAI nscodon_pdf,
                              RAI2 nsc4_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  T nuc_pdf[NUC::SIZE];
  nscodon_pdf_2_nuc_pdf(nscodon_pdf,nuc_pdf);
  nscodon_pdf_nuc_pdf_2_nsc4_pdf(nscodon_pdf,nuc_pdf,nsc4_pdf);
}




template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nsc4_pdf_pos(const RAI nscodon_pdf,
                           const unsigned pos,
                           RAI2 nsc4_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  T nuc_pdf[NUC::SIZE];
  nscodon_pdf_2_nuc_pdf_pos(nscodon_pdf,pos,nuc_pdf);
  nscodon_pdf_nuc_pdf_2_nsc4_pdf(nscodon_pdf,nuc_pdf,nsc4_pdf);
}




template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nsc4_pre_pdf(const RAI nscodon_pdf,
                           RAI2 nsc4_pdf){
  nscodon_pdf_2_nsc4_pdf_pos(nscodon_pdf,2,nsc4_pdf);
}




template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nsc4_post_pdf(const RAI nscodon_pdf,
                           RAI2 nsc4_pdf){
  nscodon_pdf_2_nsc4_pdf_pos(nscodon_pdf,0,nsc4_pdf);
}




template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nsc5_pdf(const RAI nscodon_pdf,
                       RAI2 nsc5_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  T nuc0_pdf[NUC::SIZE];
  T nuc2_pdf[NUC::SIZE];
  nscodon_pdf_2_nuc_pdf_pos(nscodon_pdf,0,nuc0_pdf);
  nscodon_pdf_2_nuc_pdf_pos(nscodon_pdf,2,nuc2_pdf);

  for(unsigned i(0);i<NSC5::SIZE;++i){
    NUC::index_t n_pre,n_post;
    NSCODON::index_t c;
    NSC5::decode(n_pre,n_post,c,static_cast<NSC5::index_t>(i));
    *(nsc5_pdf+i) = nuc2_pdf[n_pre] * *(nscodon_pdf+c) * nuc0_pdf[n_post];
  }

  bioseq_util_pdistro_check(nsc5_pdf,NSC5::SIZE);
}




template <typename FloatType,typename RAI,typename RAI2>
void
nscodon_pdf_conditioned_dinuc_pdf_2_nsc4_pdf_pos(const RAI nscodon_pdf,
                                                 const FloatType conditioned_dinuc_pdf[NUC::SIZE][NUC::SIZE],
                                                 const unsigned pos,
                                                 RAI2 nsc4_pdf){

  NUC::index_t n[CODON::BASE_SIZE];
  for(unsigned i(0);i<NSC4::SIZE;++i){
    const NSC4::index_t fi(static_cast<NSC4::index_t>(i));
    const NUC::index_t ni(NSC4::decode_nuc(fi));
    const NSCODON::index_t ci(NSC4::decode_nscodon(fi));
    NSCODON::decode(n,ci);
    *(nsc4_pdf+i) = conditioned_dinuc_pdf[n[pos]][ni] * *(nscodon_pdf+ci);
  }

  bioseq_util_pdistro_check(nsc4_pdf,NSC4::SIZE);
}




template <typename RAI,typename RAI2>
void
nscodon_pdf_dinuc_pdf_2_nsc4_pre_pdf(const RAI nscodon_pdf,
                                     const RAI dinuc_pdf,
                                     RAI2 nsc4_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  // get p(nuc_2|nuc_0)
  T n2_cond_on_n0[NUC::SIZE][NUC::SIZE];
  dinuc_pdf_2_conditioned_dinuc_pdf(dinuc_pdf,1,n2_cond_on_n0);
  nscodon_pdf_conditioned_dinuc_pdf_2_nsc4_pdf_pos(nscodon_pdf,n2_cond_on_n0,0,nsc4_pdf);
}




template <typename RAI,typename RAI2>
void
nscodon_pdf_dinuc_pdf_2_nsc4_post_pdf(const RAI nscodon_pdf,
                                      const RAI dinuc_pdf,
                                      RAI2 nsc4_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  // get p(nuc_0|nuc_2)
  T n0_cond_on_n2[NUC::SIZE][NUC::SIZE];
  dinuc_pdf_2_conditioned_dinuc_pdf(dinuc_pdf,0,n0_cond_on_n2);
  nscodon_pdf_conditioned_dinuc_pdf_2_nsc4_pdf_pos(nscodon_pdf,n0_cond_on_n2,2,nsc4_pdf);
}




template <typename RAI,typename RAI2,typename RAI3>
void
nscodon_pdf_dinuc_pdf_2_nsc5_pdf(const RAI nscodon_pdf,
                                 const RAI2 dinuc_pdf,
                                 RAI3 nsc5_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  // get p(nuc_2|nuc_0)
  T n2_cond_on_n0[NUC::SIZE][NUC::SIZE];
  dinuc_pdf_2_conditioned_dinuc_pdf(dinuc_pdf,1,n2_cond_on_n0);

  // get p(nuc_0|nuc_2)
  T n0_cond_on_n2[NUC::SIZE][NUC::SIZE];
  dinuc_pdf_2_conditioned_dinuc_pdf(dinuc_pdf,0,n0_cond_on_n2);

  NUC::index_t n[CODON::BASE_SIZE];
  for(unsigned i(0);i<NSC5::SIZE;++i){
    const NSC5::index_t fi(static_cast<NSC5::index_t>(i));
    const NUC::index_t n_pre(NSC5::decode_nuc1(fi));
    const NUC::index_t n_post(NSC5::decode_nuc2(fi));
    const NSCODON::index_t ci(NSC5::decode_nscodon(fi));
    NSCODON::decode(n,ci);
    *(nsc5_pdf+i) = n2_cond_on_n0[n[0]][n_pre] * *(nscodon_pdf+ci) * n0_cond_on_n2[n[2]][n_post];
  }

  bioseq_util_pdistro_check(nsc5_pdf,NSC5::SIZE);
}




template <typename RAI,typename RAI2>
void
nsc4_pdf_2_nuc_pdf(const RAI nsc4_pdf,
                   RAI2 nuc_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  T nscodon_pdf[NSCODON::SIZE];
  nsc4_pdf_2_nscodon_pdf(nsc4_pdf,nscodon_pdf);
  nscodon_pdf_2_nuc_pdf(nscodon_pdf,nuc_pdf);
}




template <typename RAI,typename RAI2>
void
nsc4_pdf_2_nuc_pdf_all_pos(const RAI nsc4_pdf,
                           RAI2 nuc_pdf){
  typedef typename std::iterator_traits<RAI>::value_type T;

  T nscodon_pdf[NSCODON::SIZE];
  nsc4_pdf_2_nscodon_pdf(nsc4_pdf,nscodon_pdf);
  nscodon_pdf_2_nuc_pdf_all_pos(nscodon_pdf,nuc_pdf);
}




template <typename RAI,typename RAI2>
void
nsc4_pdf_2_nscodon_pdf(const RAI nsc4_pdf,
                       RAI2 nscodon_pdf){

  for(unsigned i(0);i<NSCODON::SIZE;++i) *(nscodon_pdf+i) = 0.;

  for(unsigned i(0);i<NSC4::SIZE;++i){
    const NSCODON::index_t c1 = NSC4::decode_nscodon(static_cast<NSC4::index_t>(i));
    *(nscodon_pdf+c1) += *(nsc4_pdf+i);
  }

  bioseq_util_pdistro_check(nscodon_pdf,NSCODON::SIZE);
}



template <typename RAI,typename RAI2>
void
nsc4_pdf_2_conditional_nscodon_pdf(const RAI nsc4_pdf,
                                   const NUC::index_t nx,
                                   RAI2 nscodon_pdf){

  for(unsigned i(0);i<NSCODON::SIZE;++i) *(nscodon_pdf+i) = 0.;

  for(unsigned i(0);i<NSC4::SIZE;++i){
    const NSCODON::index_t c(NSC4::decode_nscodon(static_cast<NSC4::index_t>(i)));
    const NUC::index_t n(NSC4::decode_nuc(static_cast<NSC4::index_t>(i)));
    if(n!=nx) continue;
    *(nscodon_pdf+c) += *(nsc4_pdf+i);
  }

  pdistro_norm(nscodon_pdf,nscodon_pdf+NSCODON::SIZE);
}



template <typename RAI,typename RAI2>
void
nsc4_pdf_2_codon_boundary_dinuc_pdf(const RAI nsc4_pdf,
                                    const bool is_c4_pre,
                                    RAI2 dinuc_pdf){

  for(unsigned i(0);i<DINUC::SIZE;++i) *(dinuc_pdf+i) = 0.;

  NUC::index_t nx;
  NUC::index_t n[CODON::BASE_SIZE];
  NSCODON::index_t c;
  for(unsigned i(0);i<NSC4::SIZE;++i){
    NSC4::decode(nx,c,static_cast<NSC4::index_t>(i));
    NSCODON::decode(n,c);
    if(is_c4_pre){
      *(dinuc_pdf+DINUC::encode(nx,n[0])) += *(nsc4_pdf+i);
    } else {
      *(dinuc_pdf+DINUC::encode(n[2],nx)) += *(nsc4_pdf+i);
    }
  }

  bioseq_util_pdistro_check(dinuc_pdf,DINUC::SIZE);
}




template <typename RAI,typename RAI2>
void
nsc5_pdf_2_nscodon_pdf(const RAI nsc5_pdf,
                       RAI2 nscodon_pdf){

  for(unsigned i(0);i<NSCODON::SIZE;++i) *(nscodon_pdf+i) = 0.;

  for(unsigned i(0);i<NSC5::SIZE;++i){
    const NSCODON::index_t c1 = NSC5::decode_nscodon(static_cast<NSC5::index_t>(i));
    *(nscodon_pdf+c1) += *(nsc5_pdf+i);
  }

  bioseq_util_pdistro_check(nscodon_pdf,NSCODON::SIZE);
}




template <typename RAI,typename RAI2>
void
nsc5_pdf_2_nsc4_pdf(const RAI nsc5_pdf,
                    const bool is_c4_pre,
                    RAI2 nsc4_pdf){

  for(unsigned i(0);i<NSC4::SIZE;++i) *(nsc4_pdf+i) = 0.;

  NUC::index_t n_pre,n_post,nx;
  NSCODON::index_t c;
  for(unsigned i(0);i<NSC5::SIZE;++i){
    NSC5::decode(n_pre,n_post,c,static_cast<NSC5::index_t>(i));
    if(is_c4_pre) nx = n_pre;
    else          nx = n_post;
    *(nsc4_pdf+NSC4::encode(nx,c)) += *(nsc5_pdf+i);
  }

  bioseq_util_pdistro_check(nsc4_pdf,NSC4::SIZE);
}




template <typename RAI,typename RAI2>
void
nsc5_pdf_2_codon_boundary_dinuc_pdf(const RAI nsc5_pdf,
                                    RAI2 dinuc_pdf){

  for(unsigned i(0);i<DINUC::SIZE;++i) *(dinuc_pdf+i) = 0.;

  NUC::index_t n_pre,n_post;
  NSCODON::index_t c;
  NUC::index_t n[CODON::BASE_SIZE];
  for(unsigned i(0);i<NSC5::SIZE;++i){
    NSC5::decode(n_pre,n_post,c,static_cast<NSC5::index_t>(i));
    NSCODON::decode(n,c);
    *(dinuc_pdf+DINUC::encode(n_pre,n[0])) += *(nsc5_pdf+i);
    *(dinuc_pdf+DINUC::encode(n[2],n_post)) += *(nsc5_pdf+i);
  }
  for(unsigned i(0);i<DINUC::SIZE;++i) *(dinuc_pdf+i) *= 0.5;

  bioseq_util_pdistro_check(dinuc_pdf,DINUC::SIZE);
}




template <typename RAI,typename RAI2>
void
get_nuc_pos_distro_from_site_model_distro(const SITE_MODEL::index_t sm,
                                          const RAI smp,
                                          const unsigned pos,
                                          RAI2 p){

  const unsigned n_states(SITE_MODEL::state_size(sm));

  std::fill(p,p+NUC::SIZE,0.);

  NUC::index_t sn[SITE_MODEL::MAX_BASE_SIZE];
  for(unsigned s(0);s<n_states;++s){
    SITE_MODEL::decode_nuc(sm,s,sn);
    *(p+sn[pos]) += *(smp+s);
  }

  // if smp is normalized, then p is normalized:
  bioseq_util_pdistro_check(p,NUC::SIZE);
}




template <typename RAI,typename RAI2>
void
get_dependent_nuc_pos_distro_from_site_model_distro(const SITE_MODEL::index_t sm,
                                                    const RAI smp,
                                                    const unsigned pos,
                                                    const unsigned dep_pos,
                                                    RAI2 p){

  const unsigned n_states(SITE_MODEL::state_size(sm));

  for(unsigned i(0);i<NUC::SIZE;++i){
    std::fill(p+i*NUC::SIZE,p+(i+1)*NUC::SIZE,0.);
  }

  NUC::index_t sn[SITE_MODEL::MAX_BASE_SIZE];
  for(unsigned s(0);s<n_states;++s){
    SITE_MODEL::decode_nuc(sm,s,sn);
    *(p+sn[dep_pos]*NUC::SIZE+sn[pos]) += *(smp+s);
  }

  for(unsigned i(0);i<NUC::SIZE;++i){
    pdistro_norm(p+i*NUC::SIZE,p+(i+1)*NUC::SIZE);
  }
}




template <typename RAI,typename RAI2>
void
get_dependent_nuc_pos_distro_from_site_model_distro(const SITE_MODEL::index_t sm,
                                                    const RAI smp,
                                                    const unsigned pos,
                                                    const bool is_dep_pos,
                                                    const unsigned dep_pos,
                                                    RAI2 p){
  if(is_dep_pos){
    get_dependent_nuc_pos_distro_from_site_model_distro(sm,smp,pos,dep_pos,p);
  } else {
    get_nuc_pos_distro_from_site_model_distro(sm,smp,pos,p);
    for(unsigned i(1);i<NUC::SIZE;++i){
      std::copy(p,p+NUC::SIZE,p+i*NUC::SIZE);
    }
  }
}


template <typename RAI,typename RAI2>
void
get_nuc_distro_conditioned_from_site_model_distro(const SITE_MODEL::index_t sm,
                                                  const RAI smp,
                                                  const bool is_conditioned_5p,
                                                  RAI2 p){

  const unsigned n_states(SITE_MODEL::state_size(sm));
  const unsigned base_size(SITE_MODEL::base_size(sm));
  const bool is_coding(SITE_MODEL::is_coding(sm));
  const unsigned base_unit(is_coding ? CODON::BASE_SIZE : 1);
  const unsigned scan_size(std::min(base_unit+1,base_size));

  for(unsigned i(0);i<NUC::SIZE;++i){
    std::fill(p+i*NUC::SIZE,p+(i+1)*NUC::SIZE,0.);
  }

  if(scan_size==1){
    get_nuc_pos_distro_from_site_model_distro(sm,smp,0,p);
    for(unsigned i(1);i<NUC::SIZE;++i){
      std::copy(p,p+NUC::SIZE,p+i*NUC::SIZE);
    }
  } else {
    const unsigned shift( is_conditioned_5p ? 0 : 1 );
    for(unsigned pos(shift);pos<(scan_size-1+shift);++pos){
      NUC::index_t sn[SITE_MODEL::MAX_BASE_SIZE];
      for(unsigned s(0);s<n_states;++s){
        SITE_MODEL::decode_nuc(sm,s,sn);
        *(p+sn[pos]*NUC::SIZE+sn[pos+1-(2*shift)]) += *(smp+s);
      }
    }

    for(unsigned i(0);i<NUC::SIZE;++i){
      pdistro_norm(p+i*NUC::SIZE,p+(i+1)*NUC::SIZE);
    }
  }
}


template <typename RAI,typename RAI2>
void
get_nuc_distro_conditioned_on_5p_from_site_model_distro(const SITE_MODEL::index_t sm,
                                                        const RAI smp,
                                                        RAI2 p){
  get_nuc_distro_conditioned_from_site_model_distro(sm,smp,true,p);
}


template <typename RAI,typename RAI2>
void
get_nuc_distro_conditioned_on_3p_from_site_model_distro(const SITE_MODEL::index_t sm,
                                                        const RAI smp,
                                                        RAI2 p){
  get_nuc_distro_conditioned_from_site_model_distro(sm,smp,false,p);
}


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

// $Id: bioseq_util_pdf_conversions.h 1056 2007-12-09 22:32:55Z ctsa $

/// \file 

#ifndef __BIOSEQ_UTIL_PDF_CONVERSIONS_H
#define __BIOSEQ_UTIL_PDF_CONVERSIONS_H

#include "bioseq_util.h"
#include "../general/die.h"


namespace BIO_PDISTRO {
  enum index_t { NONE, NUC, DINUC, TRINUC, CODON, NSCODON, AA, NSAA, NSC4PRE, NSC4POST, NSC5,
                 NSC5PRE, NSC5POST,
                 NSCODON_DINUC, NUC_CODON_POS };

  inline
  index_t
  convert_from_site_model(SITE_MODEL::index_t i){
    switch(i){
    case SITE_MODEL::NUC: return NUC;
    case SITE_MODEL::DINUC: return DINUC;
    case SITE_MODEL::TRINUC: return TRINUC;
    case SITE_MODEL::CODON: return CODON;
    case SITE_MODEL::NSCODON: return NSCODON;
    case SITE_MODEL::AA: return AA;
    case SITE_MODEL::NSAA: return NSAA;
    case SITE_MODEL::NSC4PRE: return NSC4PRE;
    case SITE_MODEL::NSC4POST: return NSC4POST;
    case SITE_MODEL::NSC5: return NSC5;
    case SITE_MODEL::NSC5PRE: return NSC5PRE;
    case SITE_MODEL::NSC5POST: return NSC5POST;
    default: return NONE;
    }
  }

  inline
  SITE_MODEL::index_t
  convert_to_site_model(index_t i){
    switch(i){
    case NUC:      return SITE_MODEL::NUC;
    case DINUC:    return SITE_MODEL::DINUC;
    case TRINUC:   return SITE_MODEL::TRINUC;
    case CODON:    return SITE_MODEL::CODON;
    case NSCODON:  return SITE_MODEL::NSCODON;
    case AA:       return SITE_MODEL::AA;
    case NSAA:     return SITE_MODEL::NSAA;
    case NSC4PRE:  return SITE_MODEL::NSC4PRE;
    case NSC4POST: return SITE_MODEL::NSC4POST;
    case NSC5:     return SITE_MODEL::NSC5;
    case NSC5PRE:  return SITE_MODEL::NSC5PRE;
    case NSC5POST: return SITE_MODEL::NSC5POST;
    default: return SITE_MODEL::NONE;
    }
  }

  // dependent parameter size:
  inline
  unsigned
  param_size(const index_t i){

    const SITE_MODEL::index_t sm(convert_to_site_model(i));
    if(sm != SITE_MODEL::NONE) {
      return SITE_MODEL::state_size(sm);
    } else {
      switch(i){
      case NSCODON_DINUC: return param_size(NSCODON)+param_size(DINUC);
      case NUC_CODON_POS: return param_size(NUC)*CODON::BASE_SIZE;
      default:
        die("unknown DISTRO_TYPE");
        return 0;
      }
    }
  }
}


/// this function programatically selects most of the conversion options below
///
template <typename RAI,
          typename RAI2>
void
bio_pdistro_convert(const BIO_PDISTRO::index_t from_index,
                    const BIO_PDISTRO::index_t to_index,
                    const RAI from_pdf,
                    RAI2 to_pdf);

template <typename RAI,
          typename RAI2>
void
nuc_pdf_2_dinuc_pdf(const RAI nuc_pdf,//[NUC::SIZE],
                    RAI2 dinuc_pdf);//[DINUC::SIZE]);

template <typename RAI,typename RAI2>
void
nuc_pdf_2_codon_pdf(const RAI nuc_pdf, //[NUC::SIZE],
                    RAI2 codon_pdf); // [CODON::SIZE]);

template <typename RAI,typename RAI2>
void
nuc_pdf_2_trinuc_pdf(const RAI nuc_pdf, //[NUC::SIZE],
                     RAI2 trinuc_pdf){ //[TRINUC::SIZE]){
  nuc_pdf_2_codon_pdf(nuc_pdf,trinuc_pdf);
}

template <typename RAI,typename RAI2>
void
nuc_pdf_2_nscodon_pdf(const RAI nuc_pdf,//[NUC::SIZE],
                      RAI2 nscodon_pdf);//[NSCODON::SIZE]);

template <typename RAI,typename RAI2>
void
nuc_pdf_2_nsc4_pdf(const RAI nuc_pdf,//[NUC::SIZE],
                   RAI2 nsc4_pdf);//[NSC4::SIZE]);

template <typename RAI,typename RAI2>
void
nuc_pdf_all_pos_2_codon_pdf(const RAI nuc_pdf,//[NUC::SIZE*CODON::BASE_SIZE],
                            RAI2 codon_pdf);//[CODON::SIZE]);

template <typename RAI,typename RAI2>
void
nuc_pdf_all_pos_2_nscodon_pdf(const RAI nuc_pdf,//[NUC::SIZE*CODON::BASE_SIZE],
                              RAI2 nscodon_pdf);//[NSCODON::SIZE]);

template <typename RAI,typename RAI2>
void
nuc_pdf_all_pos_2_nsc4_pre_pdf(const RAI nuc_pdf,//[NUC::SIZE*CODON::BASE_SIZE],
                               RAI2 nsc4_pdf);//[NSC4::SIZE]);

template <typename RAI,typename RAI2>
void
nuc_pdf_all_pos_2_nsc4_post_pdf(const RAI nuc_pdf,//[NUC::SIZE*CODON::BASE_SIZE],
                                RAI2 nsc4_pdf);//[NSC4::SIZE]);

template <typename RAI,typename RAI2>
void
dinuc_pdf_2_nuc_pdf(const RAI dinuc_pdf,//[DINUC::SIZE],
                    RAI2 nuc_pdf);//[NUC::SIZE]);

template <typename RAI,typename RAI2>
void
dinuc_pdf_2_nuc_pdf_pos1(const RAI dinuc_pdf,//[DINUC::SIZE],
                         RAI2 nuc_pdf);//[NUC::SIZE]);

template <typename RAI,typename RAI2>
void
dinuc_pdf_2_conditioned_dinuc_pdf_inplace(const unsigned condition_position,
                                          RAI2 dinuc_pdf);//[DINUC::SIZE]);

template <typename RAI,typename FloatType>
void
dinuc_pdf_2_conditioned_dinuc_pdf(const RAI dinuc_pdf,//[DINUC::SIZE],
                                  const unsigned condition_position,
                                  FloatType conditioned_dinuc_pdf[NUC::SIZE][NUC::SIZE]);//[DINUC::SIZE]);

template <typename RAI,typename RAI2>
void
codon_pdf_2_nuc_pdf_pos(const RAI codon_pdf,//[CODON::SIZE],
                        const unsigned pos,
                        RAI2 nuc_pdf);//[NUC::SIZE]);

template <typename RAI,typename RAI2>
void
trinuc_pdf_2_nuc_pdf_pos(const RAI trinuc_pdf,//[TRINUC::SIZE],
                         const unsigned pos,
                         RAI2 nuc_pdf){//[NUC::SIZE]){
  codon_pdf_2_nuc_pdf_pos(trinuc_pdf,pos,nuc_pdf);
}

template <typename RAI,typename RAI2>
void
codon_pdf_2_nuc_pdf_all_pos(const RAI codon_pdf,//[CODON::SIZE],
                            RAI2 nuc_pdf);//[NUC::SIZE*CODON::BASE_SIZE]);

template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nuc_pdf(const RAI nscodon_pdf,//[NSCODON::SIZE],
                      RAI2 nuc_pdf);//[NUC::SIZE]);

template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nuc_pdf_pos(const RAI nscodon_pdf,//[NSCODON::SIZE],
                          const unsigned pos,
                          RAI2 nuc_pdf);//[NUC::SIZE]);

template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nuc_pdf_all_pos(const RAI nscodon_pdf,//[NSCODON::SIZE],
                              RAI2 nuc_pdf);//[NUC::SIZE*CODON::BASE_SIZE]);

template <typename FloatType>
void
nscodon_pdf_2_dinuc_pdfs_conditioned_on_pos(const FloatType nscodon_pdf[NSCODON::SIZE],
                                            const int pos,
                                            FloatType conditioned_nscodon_pdf[NUC::SIZE][DINUC::SIZE]);

template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nsaa_pdf(const RAI nscodon_pdf,//[NSCODON::SIZE],
                       RAI2 nsaa_pdf);//[NSAA::SIZE]);

template <typename RAI,typename RAI2,typename RAI3>
void
nscodon_pdf_nuc_pdf_2_nsc4_pdf(const RAI nscodon_pdf,//[NSCODON::SIZE],
                               const RAI2 nuc_pdf,//[NUC::SIZE],
                               RAI3 nsc4_pdf);//[NSC4::SIZE]);

template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nsc4_pdf_allnuc(const RAI nscodon_pdf,//[NSCODON::SIZE],
                              RAI2 nsc4_pdf);//[NSC4::SIZE]);

template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nsc4_pre_pdf(const RAI nscodon_pdf,//[NSCODON::SIZE],
                           RAI2 nsc4_pdf);//[NSC4::SIZE]);

template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nsc4_post_pdf(const RAI nscodon_pdf,//[NSCODON::SIZE],
                            RAI2 nsc4_pdf);//[NSC4::SIZE]);

template <typename RAI,typename RAI2>
void
nscodon_pdf_dinuc_pdf_2_nsc4_pre_pdf(const RAI nscodon_pdf,//[NSCODON::SIZE],
                                     const RAI dinuc_pdf,//[DINUC::SIZE],
                                     RAI2 nsc4_pdf);//[NSC4::SIZE]);

template <typename RAI,typename RAI2>
void
nscodon_pdf_dinuc_pdf_2_nsc4_post_pdf(const RAI nscodon_pdf,//[NSCODON::SIZE],
                                      const RAI dinuc_pdf,//[DINUC::SIZE],
                                      RAI2 nsc4_pdf);//[NSC4::SIZE]);

template <typename RAI,typename RAI2>
void
nscodon_pdf_2_nsc5_pdf(const RAI nscodon_pdf,//[NSCODON::SIZE],
                       RAI2 nsc5_pdf);//[NSC5::SIZE]);

template <typename RAI,typename RAI2,typename RAI3>
void
nscodon_pdf_dinuc_pdf_2_nsc5_pdf(const RAI nscodon_pdf,//[NSCODON::SIZE],
                                 const RAI2 dinuc_pdf,//[DINUC::SIZE],
                                 RAI3 nsc5_pdf);//[NSC5::SIZE]);

template <typename RAI,typename RAI2>
void
nsc4_pdf_2_nuc_pdf(const RAI nsc4_pdf,//[NSC4::SIZE],
                   RAI2 nuc_pdf);//[NUC::SIZE]);

template <typename RAI,typename RAI2>
void
nsc4_pdf_2_nuc_pdf_all_pos(const RAI nsc4_pdf,//[NSC4::SIZE],
                           RAI2 nuc_pdf);//[NUC::SIZE*CODON::BASE_SIZE]);

template <typename RAI,typename RAI2>
void
nsc4_pdf_2_nscodon_pdf(const RAI nsc4_pdf,//[NSC4::SIZE],
                       RAI2 nscodon_pdf);//[NSCODON::SIZE]);

template <typename RAI,typename RAI2>
void
nsc4_pdf_2_conditional_nscodon_pdf(const RAI nsc4_pdf,//[NSC4::SIZE],
                                   const NUC::index_t nx,
                                   RAI2 nscodon_pdf);//[NSCODON::SIZE]);

template <typename RAI,typename RAI2>
void
nsc4_pdf_2_codon_boundary_dinuc_pdf(const RAI nsc4_pdf,//[NSC4::SIZE],
                                    const bool is_c4_pre,
                                    RAI2 dinuc_pdf);//[DINUC::SIZE]);

template <typename RAI,typename RAI2>
void
nsc5_pdf_2_nscodon_pdf(const RAI nsc5_pdf,//[NSC5::SIZE],
                       RAI2 nscodon_pdf);//[NSCODON::SIZE]);

template <typename RAI,typename RAI2>
void
nsc5_pdf_2_nsc4_pdf(const RAI nsc5_pdf,//[NSC5::SIZE],
                    const bool is_c4_pre,
                    RAI2 nsc4_pdf);//[NSC4::SIZE]);

template <typename RAI,typename RAI2>
void
nsc5_pdf_2_codon_boundary_dinuc_pdf(const RAI nsc5_pdf,//[NSC5::SIZE],
                                    RAI2 dinuc_pdf);//[DINUC::SIZE]);



template <typename RAI,typename RAI2>
void
get_nuc_pos_distro_from_site_model_distro(const SITE_MODEL::index_t sm,
                                          const RAI smp,
                                          const unsigned pos,
                                          RAI2 p);//[NUC::SIZE]);

template <typename RAI,typename RAI2>
void
get_dependent_nuc_pos_distro_from_site_model_distro(const SITE_MODEL::index_t sm,
                                                    const RAI smp,
                                                    const unsigned pos,
                                                    const unsigned dep_pos,
                                                    RAI2 p);//[NUC::SIZE*NUC::SIZE]);

template <typename RAI,typename RAI2>
void
get_dependent_nuc_pos_distro_from_site_model_distro(const SITE_MODEL::index_t sm,
                                                    const RAI smp,
                                                    const unsigned pos,
                                                    const bool is_dep_pos,
                                                    const unsigned dep_pos,
                                                    RAI2 p);//[NUC::SIZE*NUC::SIZE]);

template <typename RAI,typename RAI2>
void
get_nuc_distro_conditioned_on_5p_from_site_model_distro(const SITE_MODEL::index_t sm,
                                                        const RAI smp,
                                                        RAI2 p);//[NUC::SIZE*NUC::SIZE]);

template <typename RAI,typename RAI2>
void
get_nuc_distro_conditioned_on_3p_from_site_model_distro(const SITE_MODEL::index_t sm,
                                                        const RAI smp,
                                                        RAI2 p);//[NUC::SIZE*NUC::SIZE]);


#include "bioseq_util_pdf_conversions.hh"

#endif

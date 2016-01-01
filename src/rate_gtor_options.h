// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// CodeAxe : phylogenetic analysis and simulation tools
//
//   http://www.phrap.org
//
//
// Copyright 2007 Christopher T Saunders (ctsa@u.washington.edu)
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

// $Id: rate_gtor_options.h 1117 2008-01-28 19:48:39Z ctsa $

/// \file

#ifndef __RATE_GTOR_OPTIONS_H
#define __RATE_GTOR_OPTIONS_H

#include "util/bio/bioseq_util.h"


/// \todo more logical reordering of enumerator -- good time to add
/// labels too
///
namespace RATE_GTOR_MODEL {
  enum index_t { NONE,
                 BINARY,
                 NUC,
                 CODON,
                 C4PRE,
                 C4POST,
                 DINUC,
                 TRINUC,
                 C5,
                 C5PRE,
                 C5POST};

  inline
  SITE_MODEL::index_t
  convert_to_site_model(const index_t i){
    switch(i){
    case BINARY: return SITE_MODEL::BINARY;
    case NUC: return SITE_MODEL::NUC;
    case CODON: return SITE_MODEL::NSCODON;
    case C4PRE: return SITE_MODEL::NSC4PRE;
    case C4POST: return SITE_MODEL::NSC4POST;
    case DINUC: return SITE_MODEL::DINUC;
    case TRINUC: return SITE_MODEL::TRINUC;
    case C5PRE: return SITE_MODEL::NSC5PRE;
    case C5: return SITE_MODEL::NSC5;
    case C5POST: return SITE_MODEL::NSC5POST;
    default: return SITE_MODEL::NONE;
    }
  }

  inline
  unsigned
  base_size(const index_t i){
    return SITE_MODEL::base_size(convert_to_site_model(i));
  }

  inline
  unsigned
  state_size(const index_t i){
    return SITE_MODEL::state_size(convert_to_site_model(i));
  }

  inline
  unsigned
  base_size_conditioned(const index_t i){
    switch(i){
    case DINUC:
    case TRINUC: return 1;
    case C4PRE:
    case C4POST:
    case C5PRE:
    case C5:
    case C5POST: return CODON::BASE_SIZE;
    default:     return SITE_MODEL::base_size(convert_to_site_model(i));
    }
  }

  inline
  unsigned
  base_overlap_size(const index_t i){
    return SITE_MODEL::base_size(convert_to_site_model(i))-base_size_conditioned(i);
  }

  /// \brief the number of nucs leading up to the repeating unit of
  /// the site
  ///
  inline
  int
  base_size_repeat_offset(const index_t i){

    switch(i){
    case C4PRE:
    case C5:
    case DINUC:  return -1;
    case C5PRE:
    case TRINUC: return -2;
    default:     return 0;
    }
  }

  inline
  unsigned
  base_overlap_left(const index_t i){
    return -base_size_repeat_offset(i);
  }

  inline
  unsigned
  base_overlap_right(const index_t i){
    return base_size(i)-base_size_conditioned(i)-base_overlap_left(i);
  }
}



struct rate_gtor_options {
};


#endif

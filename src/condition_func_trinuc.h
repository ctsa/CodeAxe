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

// $Id: condition_func_trinuc.h 1037 2007-12-03 23:58:34Z ctsa $

/// \file

#ifndef __CONDITION_FUNC_TRINUC_H
#define __CONDITION_FUNC_TRINUC_H


#include "condition_func_overlap_base.h"
#include "util/bio/bioseq_util.h"


/// \brief condition function for overlapping noncoding 3-mers:
///
struct condition_func_trinuc : public condition_func_overlap_base {
private:

  /// for each taxa, sum all possible base nuc n0 at the leaf
  /// position for each possible ancestral value (i) given the
  /// preceding two nucleotides (n)
  ///
  virtual
  void get_leaf_leg(prob_t* legp,
                    const prob_t* tprob,
                    const unsigned conditioning_state,
                    const unsigned) const {

    const unsigned CONDITIONED_NSTATES(conditioned_state_size());
    const unsigned NSTATES(state_size());

    const NUC::index_t nx = DINUC::decode_nx(static_cast<DINUC::index_t>(conditioning_state));
    const NUC::index_t n0 = DINUC::decode_n0(static_cast<DINUC::index_t>(conditioning_state));
    for(unsigned n(0);n<CONDITIONED_NSTATES;++n){
      const TRINUC::index_t tr = TRINUC::encode(nx,n0,static_cast<NUC::index_t>(n));
      const prob_t* const tprob_vec(tprob+tr*NSTATES);
      for(unsigned i(0);i<NSTATES;++i) legp[i] += tprob_vec[i];
    }
  }

  virtual
  unsigned conditioning_state_size() const { return DINUC::SIZE; }

  virtual
  unsigned conditioned_state_size() const { return NUC::SIZE; }

  virtual
  unsigned state_2_conditioning_state(unsigned i,
                                      const unsigned) const {
    // ambiguous state maps to ambiguous state:
    if(i==state_size()) return conditioning_state_size();

    NUC::index_t nuc[CODON::BASE_SIZE];
    TRINUC::decode(nuc,static_cast<TRINUC::index_t>(i));
    return DINUC::encode(nuc[0],nuc[1]);
  }
};


#endif

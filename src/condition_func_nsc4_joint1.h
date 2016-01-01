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

// $Id: condition_func_nsc4_joint1.h 1037 2007-12-03 23:58:34Z ctsa $

/// \file

#ifndef __CONDITION_FUNC_NSC4_JOINT1_H
#define __CONDITION_FUNC_NSC4_JOINT1_H


#include "condition_func_overlap_base.h"
#include "util/bio/bioseq_util.h"


/// \brief condition function for 4-mers, conditioning is first of two
/// designed for 5-mer model approximation
///
struct condition_func_nsc4_joint1 : public condition_func_overlap_base {
private:
  virtual
  void get_leaf_leg(prob_t* legp,
                    const prob_t* tprob,
                    const unsigned conditioning_state,
                    const unsigned) const {

    const unsigned NSTATES(state_size());

    const DINUC::index_t d(static_cast<DINUC::index_t>(conditioning_state));
    const NUC::index_t nx(DINUC::decode_nx(d));
    const NUC::index_t n0(DINUC::decode_n0(d));

    NUC::index_t nuc[CODON::BASE_SIZE];

    const unsigned offset(NSC4::get_nscodon_offset(nx));
    for(unsigned c(offset);c<(offset+NSCODON::SIZE);++c){
      NSCODON::index_t cod(static_cast<NSCODON::index_t>(c-offset));
      NSCODON::decode(nuc,cod);
      if(nuc[0]!=n0) continue;
      const prob_t* const tprob_vec = tprob+c*NSTATES;
      for(unsigned i(0);i<NSTATES;++i) legp[i] += tprob_vec[i];
    }
  }

  virtual
  unsigned state_size() const { return NSC4::SIZE; }

  virtual
  unsigned conditioning_state_size() const { return DINUC::SIZE; }

  virtual
  unsigned conditioned_state_size() const { return DINUC::SIZE; }

  virtual
  unsigned state_2_conditioning_state(unsigned i,
                                      const unsigned) const {
    // ambiguous state maps to ambiguous state:
    if(i==state_size()) return conditioning_state_size();

    NUC::index_t n_pre = NSC4::decode_nuc(static_cast<NSC4::index_t>(i));
    NSCODON::index_t c = NSC4::decode_nscodon(static_cast<NSC4::index_t>(i));

    NUC::index_t nuc[CODON::BASE_SIZE];
    NSCODON::decode(nuc,c);
    return DINUC::encode(n_pre,nuc[0]);
  }
};

#endif

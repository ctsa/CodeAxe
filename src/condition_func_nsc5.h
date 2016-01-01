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

// $Id: condition_func_nsc5.h 1037 2007-12-03 23:58:34Z ctsa $

/// \file

#ifndef __CONDITION_FUNC_NSC5_H
#define __CONDITION_FUNC_NSC5_H


#include "condition_func_overlap_base.h"
#include "util/bio/bioseq_util.h"


/// \brief condition function for overlapping 5-mers consisting of
/// codon with either two preceding or two following nucleotides
///
/// transform to condition on the two nucs outside of the repeating
/// codon unit
///
struct condition_func_nsc5 : public condition_func_overlap_base {
private:

  virtual
  void get_leaf_leg(prob_t* legp,
                    const prob_t* tprob,
                    const unsigned conditioning_state,
                    const unsigned) const {

    const unsigned n_conditioned_states(conditioned_state_size());
    const unsigned n_states(state_size());

    const NUC::index_t n1(DINUC::decode_nx(static_cast<DINUC::index_t>(conditioning_state)));
    const NUC::index_t n2(DINUC::decode_n0(static_cast<DINUC::index_t>(conditioning_state)));

    const unsigned offset(NSC5::get_nscodon_offset(n1,n2));
    for(unsigned c(offset);c<(offset+n_conditioned_states);++c){
      const prob_t* const tprob_vec(tprob+c*n_states);
      for(unsigned i(0);i<n_states;++i) legp[i] += tprob_vec[i];
    }
  }

  virtual
  unsigned state_size() const { return NSC5::SIZE; }

  virtual
  unsigned conditioning_state_size() const { return DINUC::SIZE; }

  virtual
  unsigned conditioned_state_size() const { return NSCODON::SIZE; }

  virtual
  unsigned state_2_conditioning_state(unsigned i,
                                      const unsigned) const {
    // ambiguous state maps to ambiguous state:
    if(i==state_size()) return conditioning_state_size();

    const NSC5::index_t ci(static_cast<NSC5::index_t>(i));
    const NUC::index_t n1(NSC5::decode_nuc1(ci));
    const NUC::index_t n2(NSC5::decode_nuc2(ci));

    return DINUC::encode(n1,n2);
  }
};


#endif

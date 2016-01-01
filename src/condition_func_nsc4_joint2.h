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

// $Id: condition_func_nsc4_joint2.h 1037 2007-12-03 23:58:34Z ctsa $

/// \file

#ifndef __CONDITION_FUNC_NSC4_JOINT2_H
#define __CONDITION_FUNC_NSC4_JOINT2_H


#include "condition_func_overlap_base.h"
#include "util/bio/bioseq_util.h"

/// \brief condition function for 4-mers, conditioning is second of two
/// designed for 5-mer model approximation
///
struct condition_func_nsc4_joint2 : public condition_func_overlap_base {
private:
  virtual
  void get_leaf_leg(prob_t* legp,
                    const prob_t* tprob,
                    const unsigned conditioning_state,
                    const unsigned) const {

    const unsigned NSTATES(state_size());

    const NSCODON::index_t c(static_cast<NSCODON::index_t>(conditioning_state));

    for(unsigned nx(0);nx<NUC::SIZE;++nx){
      const NSC4::index_t f=NSC4::encode(static_cast<NUC::index_t>(nx),c);
      const prob_t* const tprob_vec = tprob+f*NSTATES;
      for(unsigned i(0);i<NSTATES;++i) legp[i] += tprob_vec[i];
    }
  }

  virtual
  unsigned state_size() const { return NSC4::SIZE; }

  virtual
  unsigned conditioning_state_size() const { return NSCODON::SIZE; }

  virtual
  unsigned conditioned_state_size() const { return NUC::SIZE; }

  virtual
  unsigned state_2_conditioning_state(unsigned i,
                                      const unsigned) const {
    // ambiguous state maps to ambiguous state:
    if(i==state_size()) return conditioning_state_size();

    return NSC4::decode_nscodon(static_cast<NSC4::index_t>(i));
  }
};


#endif

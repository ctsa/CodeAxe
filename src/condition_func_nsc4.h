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

// $Id: condition_func_nsc4.h 1037 2007-12-03 23:58:34Z ctsa $

/// \file

#ifndef __CONDITION_FUNC_NSC4_H
#define __CONDITION_FUNC_NSC4_H


#include "condition_func_overlap_base.h"
#include "util/bio/bioseq_util.h"

#include <cmath>


/// \brief condition function for overlapping 4-mers
///
struct condition_func_nsc4 : public condition_func_overlap_base {

  typedef condition_func_overlap_base base_t;

  explicit
  condition_func_nsc4(const bool is_sqrt_result = false) : _is_sqrt_result(is_sqrt_result) {}

  virtual
  smlfloat
  transform_prob(const smlfloat p) const {
    if(_is_sqrt_result) return std::sqrt(p);
    else                return p;
  }

private:
  /// for each taxa, sum all possible codons (c) at the leaf position
  /// for each possible ancestral value (i) given a preceding nucleotide
  /// (n)
  ///
  virtual
  void get_leaf_leg(prob_t* legp,
                    const prob_t* tprob,
                    const unsigned conditioning_state,
                    const unsigned) const{

    const unsigned n_conditioned_states(conditioned_state_size());
    const unsigned n_states(state_size());

    const NUC::index_t n(static_cast<NUC::index_t>(conditioning_state));

    const unsigned offset(NSC4::get_nscodon_offset(n));
    for(unsigned c(offset);c<(offset+n_conditioned_states);++c){
      const prob_t* const tprob_vec(tprob+c*n_states);
      for(unsigned i(0);i<n_states;++i) legp[i] += tprob_vec[i];
    }
  }

  virtual
  unsigned conditioning_state_size() const { return NUC::SIZE; }

  virtual
  unsigned conditioned_state_size() const { return NSCODON::SIZE; }

  virtual
  unsigned state_2_conditioning_state(unsigned i,
                                      const unsigned) const {
    // ambiguous state maps to ambiguous state:
    if(i==state_size()) return conditioning_state_size();

    return NSC4::decode_nuc(static_cast<NSC4::index_t>(i));
  }


  bool _is_sqrt_result;
};


#endif

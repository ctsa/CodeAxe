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

// $Id: condition_func_nsc4_bidir.h 1037 2007-12-03 23:58:34Z ctsa $

/// \file

#ifndef __CONDITION_FUNC_NSC4_BIDIR_H
#define __CONDITION_FUNC_NSC4_BIDIR_H


#include "condition_func_overlap_base.h"
#include "util/bio/bioseq_util.h"
#include "util/general/die.h"

#include <cmath>


/// \brief condition function for overlapping 4-mers, using geometric
/// mean of forwards and backwards conditioning
struct condition_func_nsc4_bidir : public condition_func_overlap_base {

  typedef condition_func_overlap_base base_t;

  explicit
  condition_func_nsc4_bidir(const SITE_MODEL::index_t sm) : base_t(), _sm(sm) {
    if(_sm!=SITE_MODEL::NSC4PRE && _sm!=SITE_MODEL::NSC4POST){
      die("invalid site model for cf nsc4_bidir");
    }
  }

private:

  virtual bool is_bidirectional() const { return true; }

  ///
  virtual
  void get_leaf_leg(prob_t* legp,
                    const prob_t* tprob,
                    const unsigned conditioning_state,
                    const unsigned condition_no) const {

    const unsigned n_conditioned_states(conditioned_state_size());
    const unsigned n_states(state_size());

    const NUC::index_t n(static_cast<NUC::index_t>(conditioning_state));

    NUC::index_t nuc[SITE_MODEL::MAX_BASE_SIZE];
    NUC::index_t* nuc_ptr;

    if(condition_no==0){
      nuc[0]=n;
      nuc_ptr=nuc+1;
    } else {
      nuc[3]=n;
      nuc_ptr=nuc;
    }

    for(unsigned c(0);c<n_conditioned_states;++c){
      SITE_MODEL::decode_nuc(SITE_MODEL::TRINUC,c,nuc_ptr);
      const unsigned c4i(SITE_MODEL::encode_nuc(_sm,nuc));

      if(c4i == SITE_MODEL::ambig_state(_sm)) continue;

      const prob_t* const tprob_vec(tprob+c4i*n_states);
      for(unsigned i(0);i<n_states;++i) legp[i] += tprob_vec[i];
    }
  }

  virtual
  unsigned state_size() const { return NSC4::SIZE; }

  virtual
  unsigned conditioning_state_size() const { return NUC::SIZE; }

  virtual
  unsigned conditioned_state_size() const { return TRINUC::SIZE; }

  virtual
  unsigned state_2_conditioning_state(unsigned i,
                                      const unsigned condition_no) const {
    // ambiguous state maps to ambiguous state:
    if(i==state_size()) return conditioning_state_size();

    NUC::index_t nuc[SITE_MODEL::MAX_BASE_SIZE];
    SITE_MODEL::decode_nuc(_sm,i,nuc);
    if(condition_no==0) { return nuc[0]; }
    else                { return nuc[3]; }
  }


  const SITE_MODEL::index_t _sm;
};


#endif

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

// $Id: rate_gtor_nuc_options.h 994 2007-11-01 22:33:01Z ctsa $

/// \file

#ifndef __RATE_GTOR_NUC_OPTIONS_H
#define __RATE_GTOR_NUC_OPTIONS_H


#include <vector>


namespace RATE_MODEL_NUC {
  enum index_t { JC69,
                 K80,
                 REV,
                 NONREV,
                 F81,
                 HKY85,
                 GY94};
}

namespace CONTEXT_MODEL_NUC {
  enum index_t { INDY,
                 PRE_DOUBLET,
                 POST_DOUBLET,
                 FACTORED_TRIPLET,
                 TRIPLET,
                 CPG_ONLY,
                 TPA_ONLY,
                 CPG_2TI,
                 CPG_NONREV,
                 CPG_1TI};
}


// base class contains values directly stored in rate_gtor_nuc derived
// class contains other settings which are used only within the
// rate_gtor_nscodon_base ctor -- caller needn't be aware of this
// distinction
//
struct rate_gtor_nuc_options_base {

  rate_gtor_nuc_options_base() :
    is_use_edge_correction(true),
    _rate_model(1,RATE_MODEL_NUC::NONREV),
    _context_model(1,CONTEXT_MODEL_NUC::INDY) {}

  void set_n_site_mutation_model_cats(const unsigned i) {
    _rate_model.resize(i,_rate_model.back());
    _context_model.resize(i,_context_model.back());
  }

  RATE_MODEL_NUC::index_t rate_model(const unsigned i) const {
    return _rate_model[i];
  }

  CONTEXT_MODEL_NUC::index_t context_model(const unsigned i) const {
    return _context_model[i];
  }

  void set_rate_model(const unsigned i,
                      const RATE_MODEL_NUC::index_t rm){
    _rate_model[i] = rm;
  }

  void set_context_model(const unsigned i,
                         const CONTEXT_MODEL_NUC::index_t cm){
    _context_model[i] = cm;
  }

  void set_all_rate_model(const RATE_MODEL_NUC::index_t rm) {
    std::fill(_rate_model.begin(),_rate_model.end(),rm);
  }

  void set_all_context_model(const CONTEXT_MODEL_NUC::index_t cm) {
    std::fill(_context_model.begin(),_context_model.end(),cm);
  }


  bool is_use_edge_correction;

private:
  std::vector<RATE_MODEL_NUC::index_t> _rate_model;
  std::vector<CONTEXT_MODEL_NUC::index_t> _context_model;
};


struct rate_gtor_nuc_lock_options {
  rate_gtor_nuc_lock_options() :
    is_lock_site_rate_cats(false),
    is_lock_group_rate_cats(false) {}

  bool is_lock_site_rate_cats;
  bool is_lock_group_rate_cats;
};


/// \brief nuc mutation options passed to ctor
///
struct rate_gtor_nuc_options : public rate_gtor_nuc_options_base {
  rate_gtor_nuc_options() :
    rate_gtor_nuc_options_base(), lopt(), is_free_edge_correction(false) {}

  rate_gtor_nuc_lock_options lopt;
  bool is_free_edge_correction;
};

#endif

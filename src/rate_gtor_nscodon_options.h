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

// $Id: rate_gtor_nscodon_options.h 1146 2008-03-01 00:21:53Z ctsa $

/// \file

#ifndef __RATE_GTOR_NSCODON_OPTIONS_H
#define __RATE_GTOR_NSCODON_OPTIONS_H

#include <algorithm>
#include <vector>


// nscodon opts
//
namespace SELECT_MODEL {
  enum index_t { NONE,
                 SINGLE,
                 SYMM,
                 ASYMM,
                 FROM,
                 TO,
                 FROM_TO,
                 FROM_SYMM,
                 SYMM_TO,
                 FROM_SYMM_TO,
                 HP };
}

namespace CODON_BIAS_MODEL {
  enum index_t { NONE,
                 SYNON_RATIO,
                 SYNON_NORM_RATIO,
                 SYNON_NORM_ABS };
}


// base class contains values directly stored in rate_gtor_nscodon.
// derived class contains other settings which are used only within
// the rate_gtor_nscodon_base ctor -- caller needn't be aware of this
// distinction
//
struct rate_gtor_nscodon_options_base {

  rate_gtor_nscodon_options_base() :
    codon_bias_model(CODON_BIAS_MODEL::NONE),
    _select_model(1,SELECT_MODEL::ASYMM){}

  void set_n_select_matrix_cats(const unsigned i) {
    _select_model.resize(i,_select_model.back());
  }

  SELECT_MODEL::index_t select_model(const unsigned i) const {
    return _select_model[i];
  }

  void set_select_model(const unsigned i,
                        const SELECT_MODEL::index_t sm){
    _select_model[i] = sm;
  }

  void set_all_select_model(const SELECT_MODEL::index_t sm) {
    std::fill(_select_model.begin(),_select_model.end(),sm);
  }

  CODON_BIAS_MODEL::index_t codon_bias_model;

private:
  std::vector<SELECT_MODEL::index_t> _select_model;
};


struct rate_gtor_nscodon_lock_options {
  rate_gtor_nscodon_lock_options() :
    is_lock_site_select_cats(false),
    is_lock_group_select_cats(false) {}

  bool is_lock_site_select_cats;
  bool is_lock_group_select_cats;
};


/// \brief ctor options for rate_gtor_nscodon_base
///
struct rate_gtor_nscodon_options : public rate_gtor_nscodon_options_base {
  rate_gtor_nscodon_options() :
    rate_gtor_nscodon_options_base(), lopt() {}

  rate_gtor_nscodon_lock_options lopt;
};



// nsc5 opts
//
namespace C5_APPROX {
  enum index_t { NONE,
                 C4_JOINT,
                 C4_GMEAN };
}

#endif

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

// $Id: rate_gtor_binary.h 1109 2008-01-25 03:38:11Z ctsa $

/// \file

#ifndef __RATE_GTOR_BINARY_H
#define __RATE_GTOR_BINARY_H

#include "rate_gtor.h"
#include "util/bio/bioseq_util.h"
#include "util/math/matrix_util_io.h"
#include "util/math/prob_util_io.h"

#include <iostream>
#include <string>


/// \brief a two-state rate generator (for debugging)
///
struct rate_gtor_binary : public rate_gtor {
  typedef rate_gtor base_t;
  typedef rate_gtor_binary self_t;

  enum bogus_t { RSIZE = 2 };

  rate_gtor_binary(const rate_gtor_options& ropt,
                   const rate_gtor_sml_share& smls) : base_t(ropt,smls){
    resize(1);
  }

  rate_gtor_binary(const self_t& s,
                   const rate_gtor_sml_share& smls) : base_t(s,smls) {}

  virtual
  rate_gtor* clone(const rate_gtor_sml_share& smls) const { return new self_t(*this,smls);}

  virtual
  void
  store_state(std::ostream& os) const {
    base_t::store_state(os);
    os << "rate_gtor_binary binary_mut " << _is_train_param[0] << " " << _param[0] << "\n";
  }

  virtual
  void
  load_state(std::istream& is) {

    base_t::load_state(is);

    std::string dummy;

    bool tmpb;
    is >> dummy >> dummy >> tmpb >> _param[0];
    _is_train_param[0] = tmpb;
    fix_input_param();
  }

  virtual
  void
  rates(smlfloat rts[],
        const rates_func_options&) const {
    rts[0] = -_param[0];
    rts[1] =  _param[0];
    rts[2] =  1.;
    rts[3] = -1.;
  }
#if 0
  smlfloat
  rate_scale() const {
    return _param[0]*bg_pdistro()[0]+bg_pdistro()[1];
  }

  /// \brief report on the nucleotide model parameters:
  void
  report(std::ostream& os) const {

    const static char syms[] = "AB";

    smlfloat rts[RSIZE*RSIZE] = { 0., _param[0], 1., 0. };

    const smlfloat rate_norm(1./rate_scale());
    matrix_scale(rts,RSIZE,rate_norm);
    os <<"MUT RATE\n";
    matrix_report(rts,RSIZE,syms,os);
  }
#endif
  void
  report(std::ostream&) const {};
};


#endif

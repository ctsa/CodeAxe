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

// $Id: subs_ml_model_ci_info.h 755 2007-08-17 02:31:46Z ctsa $

/// \file

#ifndef __SUBS_ML_MODEL_CI_INFO_H
#define __SUBS_ML_MODEL_CI_INFO_H

#include "audit_info.h"
#include "subs_ml_types.h"

#include <iosfwd>
#include <string>
#include <vector>


struct param_ci_info {
  param_ci_info()
    : index(0), mean(0), variance(0), ci_plus(0), ci_minus(0), is_calculating(false) {}

  void reset() {
    index=0; mean=0; variance=0; ci_plus=0; ci_minus=0;
    points.clear();
    is_calculating=false;
  }

  void store_state(std::ostream& os) const;

  bool load_state(std::istream& is);

  unsigned index;
  smlfloat mean;
  smlfloat variance;
  smlfloat ci_plus;
  smlfloat ci_minus;
  std::vector<std::pair<smlfloat,smlfloat> > points; // stores (param_val,lhood) pairs used to calculate the ci
  bool is_calculating;
};



inline
bool operator<(const param_ci_info& lhs,
               const param_ci_info& rhs) {
  return lhs.index < rhs.index;
}




struct subs_ml_model_ci_info {

  subs_ml_model_ci_info() : _ai() {}

  explicit
  subs_ml_model_ci_info(const audit_info& ai) : _ai(ai) {}

  void store_state(std::ostream& os) const;

  void load_state(std::istream& is);

  // is stream signature correct?
  static
  bool check_state(std::istream& is);

  smlfloat alpha;
  smlfloat base_lnp; // lnp of unperturbed ml model
  std::vector<param_ci_info> val;

private:
  const audit_info _ai;
};

#endif

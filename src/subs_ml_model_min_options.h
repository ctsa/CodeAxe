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

// $Id: subs_ml_model_min_options.h 1160 2008-03-20 19:18:02Z ctsa $

/// \file

#ifndef __SUBS_ML_MODEL_MIN_OPTIONS_H
#define __SUBS_ML_MODEL_MIN_OPTIONS_H

#include "subs_ml_types.h"

#include <iosfwd>
#include <string>


namespace MIN {
  enum index_t { CONJ_DIR,
                 CONJ_GRAD,
                 CG_CD_HYBRID,
                 PRAXIS,
                 CG_PRAXIS_HYBRID,
                 EM};
}

namespace OBS_UPDATE_MODE {
  enum index_t { MAXML,
                 EQUIL};
}


extern const smlfloat DEFAULT_CONVERGE_TOLERANCE;

extern const unsigned DEFAULT_CONJ_GRAD_START_ITER;

struct subs_ml_model_min_options {
  subs_ml_model_min_options() { init(); }

  explicit
  subs_ml_model_min_options(std::istream& is);

  void store_state(std::ostream& os) const;

private:
  void init(){
    converge_value=DEFAULT_CONVERGE_TOLERANCE;
    min=MIN::CG_CD_HYBRID;
    is_prefer_diag_expm=false;
    is_cat_em=false;
    is_rootcycle=false;
    is_pingpong=false;
    max_steps=0;
    is_extra_warnings=true;
    obs_update_mode=OBS_UPDATE_MODE::MAXML;
  }

public:
  smlfloat converge_value;  // search
  MIN::index_t min;         // search
  bool is_prefer_diag_expm; // method

  // non-persistent data:
  bool is_cat_em;           // method
  bool is_rootcycle;       // method
  bool is_pingpong;         // fluff
  std::string pingpongfile; // fluff
  unsigned max_steps;       // search
  bool is_extra_warnings;    // fluff
  OBS_UPDATE_MODE::index_t obs_update_mode;
};


#endif

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

// $Id: rate_gtor.cc 1183 2008-03-27 02:12:11Z ctsa $

/// \file

#include "cat_manager.h"
#include "condition_func.h"
#include "rate_gtor.h"
#include "subs_ml_ptol.h"
#include "substk_exception.h"
#include "util/general/die.h"
#include "util/general/io_util.h"
#include "util/general/log.h"
#include "util/math/matrix_util_io.h"
#include "util/math/random_util.h"
#include "util/math/test_float.h"

#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>



/// \brief debugging pretty print of rates:
///
void
rate_gtor::
dump(const char* label,
     std::ostream& os) const {

  os << "rate parameters: " << label << "\n";
  for(unsigned i(0);i<param_size();++i){
    os << "[ " << i << " ] : " << _param[i] << "\n";
  }
  os << "\n";
}


void
rate_gtor::
check_param() const {
  const unsigned ps(param_size());
  for(unsigned i(0);i<ps;++i) {
    if(is_float_invalid(_param[i])){
      std::ostringstream oss;
      oss << "rate_gtor.check_param(): Invalid model parameter: " << i << " " << _param[i] ;
      throw substk_exception(oss.str().c_str());
    }
  }
}



void
rate_gtor::
reset_param_internal(const PARAM_INIT_TYPE::index_t pinit){

  // randomization parameters
  //
  static const smlfloat rate_min(0.01);
  static const smlfloat rate_max(0.66);

  // derived parameters:
  //
  static const smlfloat rate_range(rate_max-rate_min);

  if(pinit == PARAM_INIT_TYPE::RANDOM){
    const unsigned ps(param_size());
    for(unsigned i(0);i<ps;++i){
      if(_is_train_param[i]) { _param[i] = random_uniform()*rate_range+rate_min; }
      else                   { _param[i] = 1.; }
    }

  } else if(pinit == PARAM_INIT_TYPE::START){
    set_uniform_param_state(1.);

  } else { die("unknown param_init_type"); }
}



void
rate_gtor::
load_state(std::istream& is) {
  load_state_internal(is);
  _is_param_init=true;
  bg_pdistro_update_internal();
  fix_input_param();
}



void
rate_gtor::
fix_input_param() {
  if(_is_param_init){
    fix_input_param_internal();
    check_param();
  }
}



condition_func*
rate_gtor::
condition_func_factory() const {
  return new condition_func;
}

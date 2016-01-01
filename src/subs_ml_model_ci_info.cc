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

// $Id: subs_ml_model_ci_info.cc 755 2007-08-17 02:31:46Z ctsa $
//

/// \file

#include "subs_ml_model_ci_info.h"
#include "subs_ml_print_util.h"
#include "util/general/io_util.h"
#include "util/general/log.h"
#include "util/general/metatags.h"

#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <sstream>



const unsigned MODEL_CI_FILE_VERSION(3);


const char * const section_id = "subs_ml_model_ci_info";

const char * const version_id = "model_ci_file_version";


static
void
load_die(std::istream& is,
         const char* const sid,
         const std::string& tmpstr)  NORETURN_TAG;


static
void
load_die(std::istream& is,
         const char* const sid,
         const std::string& tmpstr){

  log_os << "ERROR:: -- Invalid confidence interval file format. --\n"
         << "ERROR:: expected_section: " << sid << "\n"
         << "ERROR:: section: " << tmpstr << "\n"
         << "\n";
  log_os << is.rdbuf();
  abort();
}


static
void
load_advance(std::istream& is,
             const char* const second = 0){

  std::string tmpstr;
  is >> tmpstr;
  if(tmpstr != section_id) load_die(is,section_id,tmpstr);
  is >> tmpstr;
  if(second && tmpstr != second) load_die(is,second,tmpstr);
}







void
param_ci_info::
store_state(std::ostream& os) const {

  os << std::setprecision(SUBS_ML_PRINT_MODEL_PRECISION);
  os << section_id << " freeparam " << index;
  if(! is_calculating){
    os << " " << mean << " " << variance
       << " " << ci_plus << " " << ci_minus;
  }
  const unsigned ps(points.size());
  for(unsigned i(0);i<ps;++i){
    os << " " << std::setprecision(SUBS_ML_PRINT_MODEL_PRECISION)
       << points[i].first
       << " " << std::setprecision(SUBS_ML_PRINT_LNP_PRECISION)
       << points[i].second;
  }
  os << "\n";
}





bool
param_ci_info::
load_state(std::istream& is){

  reset();
  bool b;
  std::string buf;
  if(!is) return static_cast<bool>(is);
  if((b=getline(is,buf))){
    std::istringstream sbuf(buf);
    if(!sbuf) return static_cast<bool>(sbuf);

    load_advance(sbuf,"freeparam");
    sbuf >> index;

    is_calculating=true;
    if(sbuf >> mean >> variance){
      is_calculating=false;
    }

    sbuf >> ci_plus >> ci_minus;

    smlfloat param,lnp;
    while(sbuf >> param >> lnp){
      points.push_back(std::make_pair(param,lnp));
    }
  }
  return b;
}




void
subs_ml_model_ci_info::
store_state(std::ostream& os) const {

  os << section_id << " " << version_id << " " << MODEL_CI_FILE_VERSION << "\n";
  _ai.store_state(os);

  os << std::setprecision(SUBS_ML_PRINT_MODEL_PRECISION);
  os << section_id << " alpha " << alpha << "\n";
  os << std::setprecision(SUBS_ML_PRINT_LNP_PRECISION);
  os << section_id << " base_lnp " << base_lnp << "\n";

  const unsigned s(val.size());
  for(unsigned i(0);i<s;++i){
    val[i].store_state(os);
  }
}




bool
subs_ml_model_ci_info::
check_state(std::istream& is) {

  std::string tmpstr;
  if(! (is >> tmpstr)) return false;
  if(tmpstr != section_id) load_die(is,section_id,tmpstr);
  if(! (is >> tmpstr)) return false;
  return (tmpstr == version_id);
}




void
subs_ml_model_ci_info::
load_state(std::istream& is) {

  bool is_v1_read(false);

  {  // check file version
    unsigned tmpu;
    load_advance(is,version_id);
    is >> tmpu;
    if(tmpu != MODEL_CI_FILE_VERSION){
      if(tmpu == 1){
        is_v1_read=true;
      } else {
        log_os << "ERROR:: bad model_ci_info file version number. expected: "
               << MODEL_CI_FILE_VERSION << "\n";
        exit(EXIT_FAILURE);
      }
    }
  }

  _ai.skip_state(is);

  load_advance(is,"alpha");
  is >> alpha;

  if(! is_v1_read){
    load_advance(is,"base_lnp");
    is >> base_lnp;
  } else {
    base_lnp=1; // intentionally set to an invalid value;
  }

  // clear out the end of a previous line:
  clear_whitespace_line(is);
  val.clear();

  param_ci_info pc;
  while(pc.load_state(is)){ val.push_back(pc); }
}

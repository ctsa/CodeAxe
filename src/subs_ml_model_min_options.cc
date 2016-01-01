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

// $Id: subs_ml_model_min_options.cc 1153 2008-03-18 00:24:07Z ctsa $

/// \file

#include "subs_ml_model_min_options.h"
#include "util/general/io_util.h"
#include "util/general/log.h"

#include <cassert>
#include <cstdlib>

#include <iostream>
#include <sstream>


const smlfloat DEFAULT_CONVERGE_TOLERANCE(1e-4);

const unsigned DEFAULT_CONJ_GRAD_START_ITER(10);


const char* const section_id = "subs_ml_model_min_options";

const char* const converge_label = "converge_tolerance";
const char* const min_label = "minimizer_type";
const char* const expm_label = "is_prefer_diag_expm";
const char* const obs_up_label = "obs_update_mode";

const char* const end_label = "END";



void
subs_ml_model_min_options::
store_state(std::ostream& os) const {

  os << section_id << " " << converge_label << " " << converge_value << "\n";
  os << section_id << " " << min_label << " " << min << "\n";
  os << section_id << " " << expm_label << " " << is_prefer_diag_expm << "\n";
  os << section_id << " " << obs_up_label << " " << obs_update_mode << "\n";
  os << section_id << " " << end_label << "\n";
}



subs_ml_model_min_options::
subs_ml_model_min_options(std::istream& is){

  init();

  // clear out the end of a previous line:
  clear_whitespace_line(is);

  std::string buf;

  while(getline(is,buf)){
    std::istringstream sbuf(buf);

    std::string tmpstr;
    sbuf >> tmpstr;
    if(tmpstr != section_id){
      log_os << "ERROR:: -- Invalid model file format. --\n"
             << "ERROR:: expected_section: " << section_id << "\n"
             << "ERROR:: section: " << tmpstr << "\n"
             << "ERROR:: line buffer: " << buf << "\n";
      log_os << is.rdbuf();
      abort();
    }

    sbuf >> tmpstr;

    if(tmpstr == end_label) return;

    if       (tmpstr == converge_label){
      sbuf >> converge_value;
    } else if(tmpstr == min_label){
      min = static_cast<MIN::index_t>(read_int(sbuf));
    } else if(tmpstr == expm_label){
      sbuf >> is_prefer_diag_expm;
    } else if(tmpstr == obs_up_label){
      obs_update_mode = static_cast<OBS_UPDATE_MODE::index_t>(read_int(sbuf));
    } else {
      log_os << "ERROR:: invalid model file format. subsection: " << tmpstr << "\n";
      log_os << is.rdbuf();
      abort();
    }
  }

  log_os << "ERROR:: invalid model file format. no end to section: " << section_id << "\n";
  log_os << is.rdbuf();
  abort();
}

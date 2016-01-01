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

// $Id: param_util.cc 791 2007-09-06 21:55:22Z ctsa $

/// \file

#include "param_util.h"

void
independent_to_dependent_param_block(const std::vector<bool>& is_train_param,
                                     const std::vector<bool>& is_free_block_start,
                                     const unsigned i,
                                     std::vector<smlfloat>& param){

  smlfloat block_scale(0);
  unsigned block_size(0);  // number of free dependent parameters
  unsigned block_start(0);

  for(unsigned j(i);;--j){
    if(is_train_param[j]){
      block_scale += param[j];
      block_size++;
    }
    if(is_free_block_start[j]) {
      block_start=j;
      break;
    }
    if(j==0) die("Invalid independent param block structure.");
  }

  if(block_size==0) return;

  const smlfloat A((2.*static_cast<smlfloat>(block_size))/(block_scale));

  bool is_first(true);
  for(unsigned j(block_start);j<=i;++j){
    if(is_train_param[j]){
      if(is_first){
        param[j] = A-1.;
        is_first = false;
      } else {
        param[j] = (param[j]*A)-1.;
      }

      param[j]=std::fabs(param[j]);
    }
  }
}

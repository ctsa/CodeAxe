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

// $Id: matrix_util.hh 743 2007-08-14 15:47:12Z ctsa $

/// \file

template <typename FloatType>
void
matrix_state_reduction(FloatType* reduced_mat,
                       const unsigned reduced_size,
                       const FloatType* full_mat,
                       const unsigned full_size,
                       const unsigned* reduction_map,
                       const bool is_average,
                       const bool is_exclude_diag){

  array_zero(reduced_mat,reduced_size*reduced_size);

  unsigned* reduced_add_count(0);
  if(is_average){
    reduced_add_count = new unsigned[reduced_size*reduced_size];
    array_zero(reduced_add_count,reduced_size*reduced_size);
  }

  for(unsigned i(0);i<full_size;++i){
    const unsigned ri = reduction_map[i];
    for(unsigned j(0);j<full_size;++j){
      const unsigned rj = reduction_map[j];

      if(is_exclude_diag && i==j) continue;

      reduced_mat[rj+ri*reduced_size] += full_mat[j+i*full_size];
      if(is_average) reduced_add_count[rj+ri*reduced_size]++;
    }
  }

  if(is_average){
    for(unsigned ri(0);ri<reduced_size;++ri){
      for(unsigned rj(0);rj<reduced_size;++rj){
        const unsigned val = reduced_add_count[rj+ri*reduced_size];
        if(val) reduced_mat[rj+ri*reduced_size] /= static_cast<FloatType>(val);
      }
    }
    delete [] reduced_add_count;
  }
}

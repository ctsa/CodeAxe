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

// $Id: matrix_util_io.hh 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include <algorithm>
#include <ostream>
#include <iomanip>

template <typename FloatType>
void
matrix_report(const FloatType mat[],
              const unsigned N,
              const char* syms,
              std::ostream& os,
              const unsigned cell_width,
              const unsigned prec){

  // don't try to print something too big
  const unsigned Nprint(std::min(N,static_cast<unsigned>(MAX_MATRIX_PRINT)));

  // top row labels:
  os << "  ";
  for(unsigned i(0);i<Nprint;++i) {
    for(unsigned j(0);j<cell_width;++j){
      if(j==(cell_width-3)){
        os << syms[i];
      } else {
        os << ' ';
      }
    }
  }
  os << '\n';

  for(unsigned i(0);i<Nprint;++i){
    os << syms[i] << " ";
    for(unsigned j(0);j<Nprint;++j){
      os << " ";
      FloatType val = mat[j+i*N];
      if ( val < 1.e-99 && val >= 0.){
        if( i == j ){
          for(unsigned k(0);k<(cell_width-2);++k){
            os << '-';
          }
          os << syms[i];
        } else {
          for(unsigned k(0);k<(cell_width-2);++k){
            os << ' ';
          }
          os << "0";
        }
      } else {
        os.unsetf(std::ios::fixed);
        os << std::right
           << std::setw(cell_width-1)
           << std::setprecision(prec);
        os << val;
      }
    }
    os << " " << syms[i] << "\n";
  }
  os << "\n";
}




template <typename FloatType>
void
matrix_dump(const FloatType mat[],
            const unsigned N,
            std::ostream& os,
            int prec){

  // don't try to print something too big
  unsigned Nprint = N;
  if(Nprint>MAX_MATRIX_PRINT) Nprint = MAX_MATRIX_PRINT;

  for(unsigned i(0);i<Nprint;++i){
    for(unsigned j(0);j<Nprint;++j){
      os << " ";
      FloatType val = mat[j+i*N];
//       if ( val < FloatType(1.e-99) && val > FloatType(0.)){
//         if( i == j ){
//           os << "-------";
//         } else {
//           os << "      0";
//         }
//       } else {
        os
           << std::setw(7)
           << std::setprecision(prec)
           << val;
//       }
    }
    os << "\n";
  }
  os << "\n";
}

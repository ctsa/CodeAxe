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

// $Id: prob_util_io.hh 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "prob_util.h"
#include "matrix_util.h"

#include <ostream>
#include <iomanip>


/// \brief simple debugging output
///
template <typename FloatType>
void
pdistro_dump(const FloatType pn[],
             const unsigned N,
             std::ostream& os,
             const char* label) {

  static const char foostr[] = "pdistro";

  if(!label) label = foostr;

  for(unsigned i=0;i<N;++i){
    os << label << "[" << std::setw(2) << i << "] : "
       << std::setprecision(6) << pn[i] << "\n";
  }
}


/// \brief simple debugging output
///
template <typename FloatType,
          typename SymType>
void
pdistro_report(const FloatType* pdf,
               const unsigned N,
               const SymType* syms,
               std::ostream& os){
  os << std::setprecision(5) << std::left << std::fixed;
  for(unsigned i(0);i<N;++i){
    os << syms[i] << " : " << std::setw(8) << pdf[i] << "\n";
  }
  os << "\n";
}

#if 0
/// report a transition probability matrix
///
template <typename FloatType>
void
report_transition_prob(const FloatType* tprob, // N*N
                       const char* label,
                       const char* syms,
                       std::ostream& os,
                       const unsigned N){
  // dump basic rates
  std::string s(label);
  std::string s1 = s+" : TPROB";
  os << s1 << "\n";
  matrix_report(tprob,N,syms,os);

  // dump rate asymmetries:
  std::string s2 = s+" : LLR-BITS-ASYM-TPROB log2(ij/ji)";
  FloatType* asym = new FloatType[N*N];
  matrix_l2lr_asym(asym,tprob,N);
  os << s2 << "\n";
  matrix_report(asym,N,syms,os);
  delete [] asym;
}


template <typename FloatType>
void
report_transition_prob_flux(const FloatType* tprob,        // N*N
                            const FloatType* from_distro,  // N
                            const char* label,
                            const char* syms,
                            std::ostream& os,
                            const unsigned N){
  // dump flux:
  std::string s(label);
  std::string s1 = s+" : FLUX";
  FloatType* flux = new FloatType[N*N];
  tprob_to_joint(flux,tprob,from_distro,N);
  os << s1 << "\n";
  matrix_report(flux,N,syms,os);

  // dump flux asymmetries:
  std::string s2 = s+" : LLR-BITS-ASYM-FLUX log2(ij/ji)";
  FloatType* asym = new FloatType[N*N];
  matrix_l2lr_asym(asym,flux,N);
  os << s2 << "\n";
  matrix_report(asym,N,syms,os);

  delete [] flux;
  delete [] asym;
}
#endif

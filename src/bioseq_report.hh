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

// $Id: bioseq_report.hh 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "util/bio/bioseq_util.h"

#include <cmath>

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <limits>
#include <vector>



template <typename T>
struct smc1{
  smc1(unsigned _i,
       unsigned _j,
       T _val) : i(_i),j(_j),val(_val) {}

  bool
  operator<(const smc1& right) const {
    return val < right.val;
  }

  unsigned i;
  unsigned j;
  T val;
};


template <typename FloatType>
void
report_sorted_nsaa_pair_values(const FloatType* nsaa_val,
                               std::ostream& os,
                               const char* label){

  //// dump sorted aa values:
  // put non-zero,non-diag values into a sorting struct:
  std::vector<smc1<FloatType> > vsm;

  for(unsigned i(0);i<NSAA::SIZE;++i){
    for(unsigned j(0);j<NSAA::SIZE;++j){
      if(i == j) continue;
      FloatType val = nsaa_val[j+i*NSAA::SIZE];
      if(std::fabs(val) < std::numeric_limits<smlfloat>::epsilon() ) continue;
      vsm.push_back(smc1<FloatType>(i,j,val));
    }
  }

  std::sort(vsm.begin(),vsm.end());
  std::reverse(vsm.begin(),vsm.end());

  os << "SORTED " << label << " PARAMETERS:\n";
  os << "(from->to=param)\n";
  for(unsigned i(0);i<vsm.size();++i){
    os << NSAA::syms[vsm[i].i] << "->"
       << NSAA::syms[vsm[i].j] << "="
       << std::setprecision(7)
       << std::left
       << vsm[i].val;

    { // add codon transition list:
      std::vector<nscodon_pair> vc;
      get_codon_paths(vc,
                      static_cast<NSAA::index_t>(vsm[i].i),
                      static_cast<NSAA::index_t>(vsm[i].j));

      os << " &   "
         << NSAA::syms[vsm[i].i] << "->"
         << NSAA::syms[vsm[i].j] << " codons: ";
      for(unsigned j(0);j<vc.size();++j){
        if(j) os << ",";
        os << NSCODON::print(vc[j].c1) << "->" << NSCODON::print(vc[j].c2);
      }
    }
    os << "\n";
  }
  os << "\n";
}

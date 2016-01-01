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

// $Id: rate_gtor_nscodon_base_util.cc 1145 2008-02-28 20:07:28Z ctsa $

/// \file

#include "rate_gtor_nscodon_base_util.h"



smlfloat
average_codon_selection(const rate_gtor_nscodon_base& r,
                        const NUC::index_t* nuc,
                        const NUC::index_t mut_nuc,
                        const unsigned mut_pos,
                        const rates_func_options& opt,
                        const prob_t* nscodon_pdistro){

  prob_t total_codon_prob(0.);
  smlfloat total_codon_selection(0.);
  for(unsigned c1(0); c1<NSCODON::SIZE; ++c1){
    const NSCODON::index_t ci1(static_cast<NSCODON::index_t>(c1));
    NUC::index_t cin1[CODON::BASE_SIZE];
    NSCODON::decode(cin1,ci1);

    bool is_skip(false);
    for(unsigned i(0);i<CODON::BASE_SIZE;++i){
      if(nuc[i] != NUC::N && nuc[i] != cin1[i]){
        is_skip=true;
        break;
      }
    }
    if(is_skip) continue;

    const prob_t codon_prob(nscodon_pdistro[c1]);
    total_codon_prob += codon_prob;

    NUC::index_t cin2[CODON::BASE_SIZE];
    for(unsigned i(0);i<CODON::BASE_SIZE;++i) cin2[i] = cin1[i];
    cin2[mut_pos] = mut_nuc;

    const NSCODON::index_t ci2(NSCODON::encode(cin2));

    smlfloat codon_selection(1.);
    if( ci2 == NSCODON::NNN ) { codon_selection = 0.; }// stop site
    else {

      if(! opt.is_skip_nonaa_selection){
        codon_selection *= r.codon_bias(ci1,ci2);
      }

      if(! opt.is_skip_aa_selection){
        const NSAA::index_t ci1aa(codon_trans_known(ci1));
        const NSAA::index_t ci2aa(codon_trans_known(ci2));

        if(ci1aa != ci2aa){
          codon_selection *= r.cat_nsaa_param(ci1aa,ci2aa,opt.cat,opt.branch_cat_set);
        }
      }
    }

#if 0
    log_os << NSCODON::print(s1c) << "->" << NSCODON::print(s2c)
           <<"::select: " << codon_selection << " " << codon_prob << "\n";
#endif
    total_codon_selection += codon_selection*codon_prob;
  }
  return total_codon_selection/total_codon_prob;
}

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

// $Id: tang_ei.cc 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "site_data_stats.h"
#include "tang_ei.h"
#include "util/bio/bioseq_util.h"
#include "util/general/die.h"

#include <algorithm>
#include <iostream>

void
get_ei(const site_data& sd,
       const std::string& tax_label1,
       const std::string& tax_label2){

  static const smlfloat kappa(3.4);
  static const smlfloat c_factor(0.930205);

  if(sd.sm != SITE_MODEL::NSCODON) pass_away("Evolutionary index requires NSCODON site type");

  const unsigned seqid1(sd.taxid.getid(tax_label1));
  const unsigned seqid2(sd.taxid.getid(tax_label2));

  static const unsigned nc2(NSCODON::SIZE*NSCODON::SIZE);
  static const unsigned naa2(NSAA::SIZE*NSAA::SIZE);

  smlfloat aa_obs_count[naa2];
  unsigned codon_obs_count_1d[NSCODON::SIZE];

  { // get counts of observed aa substitutions:
    // first get counts of codon substitutions:
    unsigned codon_obs_count[nc2];
    std::fill(codon_obs_count,codon_obs_count+nc2,0);

    {
      site_count_map::const_iterator i,i_end=sd.site_count.end();

      // first pass: only count codons pairs with 0 or 1 nuc differences
      for(i=sd.site_count.begin();i!=i_end;++i){
        const unsigned state1(i->first.get_taxid(seqid1));
        const unsigned state2(i->first.get_taxid(seqid2));

        if(state1>=NSCODON::SIZE) continue;
        if(state2>=NSCODON::SIZE) continue;

        const NSCODON::index_t c1(static_cast<NSCODON::index_t>(state1));
        const NSCODON::index_t c2(static_cast<NSCODON::index_t>(state2));

        NUC::index_t n1[CODON::BASE_SIZE];
        NUC::index_t n2[CODON::BASE_SIZE];

        NSCODON::decode(n1,c1);
        NSCODON::decode(n2,c2);

        unsigned n_diff(0);
        for(unsigned k(0);k<CODON::BASE_SIZE;++k) if(n1[k]!=n2[k]) n_diff++;

        if(n_diff>=2) continue;

        const unsigned cindex(c2+c1*NSCODON::SIZE);
        site_count_data_type::const_iterator j,j_end=i->second.end();

        for(j=i->second.begin();j!=j_end;++j) codon_obs_count[cindex] += j->second;
      }
#if 0
      // second pass: only count codons pairs with 2 nuc differences
      for(i=sd.site_count.begin();i!=i_end;++i){
        const unsigned state1(i->first.get_taxid(seqid1));
        const unsigned state2(i->first.get_taxid(seqid2));

        if(state1>=NSCODON::SIZE) continue;
        if(state2>=NSCODON::SIZE) continue;

        const NSCODON::index_t c1(static_cast<NSCODON::index_t>(state1));
        const NSCODON::index_t c2(static_cast<NSCODON::index_t>(state2));

        NUC::index_t n1[CODON::BASE_SIZE];
        NUC::index_t n2[CODON::BASE_SIZE];

        NSCODON::decode(n1,c1);
        NSCODON::decode(n2,c2);

        unsigned n_diff(0);
        for(unsigned k(0);k<CODON::BASE_SIZE;++k) if(n1[k]!=n2[k]) n_diff++;

        if(n_diff!=2) continue;

        unsigned d1,d2;


        const unsigned cindex(c2+c1*NSCODON::SIZE);
        site_count_data_type::const_iterator j,j_end=i->second.end();

        for(j=i->second.begin();j!=j_end;++j) codon_obs_count[cindex] += j->second;
      }
#endif
    }

    // get 1dcodon counts:
    std::fill(codon_obs_count_1d,codon_obs_count_1d+NSCODON::SIZE,0);
    for(unsigned c1(0);c1<NSCODON::SIZE;++c1){
      for(unsigned c2(0);c2<NSCODON::SIZE;++c2){
        const unsigned cindex(c2+c1*NSCODON::SIZE);
        codon_obs_count_1d[c1] += codon_obs_count[cindex];
        codon_obs_count_1d[c2] += codon_obs_count[cindex];
      }
    }

    // now reduce to aa counts and fold:
    std::fill(aa_obs_count,aa_obs_count+naa2,0);

    for(unsigned c1(0);c1<NSCODON::SIZE;++c1){
      const NSAA::index_t aa1(codon_trans_known(static_cast<NSCODON::index_t>(c1)));

      for(unsigned c2(0);c2<NSCODON::SIZE;++c2){
        const NSAA::index_t aa2(codon_trans_known(static_cast<NSCODON::index_t>(c2)));

        const unsigned cindex(c2+c1*NSCODON::SIZE);
        const unsigned aaindex(aa2+aa1*NSAA::SIZE);

        aa_obs_count[aaindex] += codon_obs_count[cindex];
      }
    }

    for(unsigned aa1(0);aa1<NSAA::SIZE;++aa1){
      for(unsigned aa2(aa1+1);aa2<NSAA::SIZE;++aa2){
        const unsigned aaindex1(aa2+aa1*NSAA::SIZE);
        const unsigned aaindex2(aa1+aa2*NSAA::SIZE);
        const smlfloat val((aa_obs_count[aaindex1]+aa_obs_count[aaindex2])/2.);
        aa_obs_count[aaindex1]=val;
        aa_obs_count[aaindex2]=val;
      }
    }
  }



  smlfloat aa_expect_count[naa2];
  { // get counts of expected aa substitutions:
    unsigned ti_codon_count[nc2];
    unsigned tv_codon_count[nc2];

    std::fill(ti_codon_count,ti_codon_count+nc2,0);
    std::fill(tv_codon_count,tv_codon_count+nc2,0);

    for(unsigned c1(0);c1<NSCODON::SIZE;++c1){
      NUC::index_t n1[CODON::BASE_SIZE];
      NSCODON::decode(n1,static_cast<NSCODON::index_t>(c1));
      for(unsigned j(0);j<CODON::BASE_SIZE;++j){
        for(unsigned k(0);k<NUC::SIZE;++k){
          const NUC::index_t nk2(static_cast<NUC::index_t>(k));
          if(nk2==n1[j]) continue;

          NUC::index_t n2[CODON::BASE_SIZE];
          for(unsigned i(0);i<CODON::BASE_SIZE;++i) n2[i]=n1[i];
          n2[j] = nk2;

          // test for stop codon:
          const NSCODON::index_t c2(NSCODON::encode(n2));
          if(c2 == NSCODON::NNN) continue;

          const unsigned cindex(c2+c1*NSCODON::SIZE);
          if(NUC::is_transition(n1[j],n2[j])) ti_codon_count[cindex] += 1;
          else      tv_codon_count[cindex] += 1;
        }
      }
    }

    // scale by known codon counts:
    for(unsigned c1(0);c1<NSCODON::SIZE;++c1){
      const unsigned scale(codon_obs_count_1d[c1]);
      for(unsigned c2(0);c2<NSCODON::SIZE;++c2){
        const unsigned cindex(c2+c1*NSCODON::SIZE);
        ti_codon_count[cindex] *= scale;
        tv_codon_count[cindex] *= scale;
      }
    }

    // now reduce to aa counts:
    unsigned ti_aa_count[naa2];
    unsigned tv_aa_count[naa2];

    std::fill(ti_aa_count,ti_aa_count+naa2,0);
    std::fill(tv_aa_count,ti_aa_count+naa2,0);

    for(unsigned c1(0);c1<NSCODON::SIZE;++c1){
      const NSAA::index_t aa1(codon_trans_known(static_cast<NSCODON::index_t>(c1)));

      for(unsigned c2(0);c2<NSCODON::SIZE;++c2){
        const NSAA::index_t aa2(codon_trans_known(static_cast<NSCODON::index_t>(c2)));

        const unsigned cindex(c2+c1*NSCODON::SIZE);
        const unsigned aaindex(aa2+aa1*NSAA::SIZE);

        ti_aa_count[aaindex] += ti_codon_count[cindex];
        tv_aa_count[aaindex] += tv_codon_count[cindex];
      }
    }

    // make final expectation matrix
    // 1. combine ti and tv
    // 2. scale to observed synonymous level
    // 3. fold off-diagonal elements for symmetry
    //
    {
      for(unsigned aa1(0);aa1<NSAA::SIZE;++aa1){
        for(unsigned aa2(aa1);aa2<NSAA::SIZE;++aa2){
          if(aa1==aa2){
            const unsigned aaindex(aa2+aa1*NSAA::SIZE);
            aa_expect_count[aaindex] = ti_aa_count[aaindex]*kappa+tv_aa_count[aaindex];
          } else {
            const unsigned aaindex1(aa2+aa1*NSAA::SIZE);
            const unsigned aaindex2(aa1+aa2*NSAA::SIZE);
            const smlfloat val((ti_aa_count[aaindex1]+ti_aa_count[aaindex2])*kappa
                               +(tv_aa_count[aaindex1]+tv_aa_count[aaindex2]));
            aa_expect_count[aaindex1] = val/2.;
            aa_expect_count[aaindex2] = val/2.;
          }
        }
      }

      // get no of expected and observed changes in alignment:
      smlfloat n_synon_exp(0);
      for(unsigned aa(0);aa<NSAA::SIZE;++aa) n_synon_exp += aa_expect_count[aa*(1+NSAA::SIZE)];

      smlfloat n_synon_obs(0);
      for(unsigned i(0);i<NSAA::SIZE;++i) n_synon_obs += aa_obs_count[i*(1+NSAA::SIZE)];

      const smlfloat scale(n_synon_obs/n_synon_exp);
      for(unsigned i(0);i<naa2;++i) aa_expect_count[i] *= scale;
    }
  }


  // report the whole shebang:
  std::ostream& os(std::cout);
  os << "aa1 aa2 obs exp c*obs/exp\n";
  for(unsigned i(0);i<NSAA::SIZE;++i){
    for(unsigned j(0);j<NSAA::SIZE;++j){
      if(i==j) continue;
      const unsigned aaindex(j+i*NSAA::SIZE);
      if(aa_expect_count[aaindex] <= 0.) continue;
      os << NSAA::syms[i] << " " << NSAA::syms[j] << " "
         << aa_obs_count[aaindex] << " " << aa_expect_count[aaindex] << " "
         << c_factor*aa_obs_count[aaindex]/aa_expect_count[aaindex] << "\n";
    }
  }
}

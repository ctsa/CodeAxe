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

// $Id: sim_context_discrete.cc 1130 2008-01-29 03:41:42Z ctsa $

/// \file

#include "sim_context_discrete.h"
#include "sim_context_shared.h"
#include "simple_util.h"
#include "util/bio/bioseq_util.h"
#include "util/general/log.h"
#include "util/math/random_util.h"

#include <ostream>



/// \todo pass n_cats in a real parameter, not this indirect cat_sim_time business...
///


static
void
simulate_steps_nuc(const std::vector<unsigned>& cat_seq,
                   const std::vector<unsigned>& group_seq,
                   std::vector<unsigned>& leaf,
                   const prob_t* const * nucsub,
                   const prob_t* const * nucsub_cdf,
                   const prob_t* const * bg_nuc_cdf_cat,
                   const unsigned nsteps,
                   const smlfloat last_fracstep,
                   std::vector<smlfloat>& cat_sim_time){

  // ring buffer (better described as a tiny fifo?) holds the +/-1
  // codon view around the current codon
  static const unsigned RING_SIZE(3);

  const unsigned n_cats(cat_sim_time.size());
  for(unsigned i(0);i<n_cats;++i) cat_sim_time[i] = 0.;

  NUC::index_t n_ring[RING_SIZE];

  const unsigned lsize(leaf.size());

  smlfloat fracstep(1.);
  for(unsigned step(0);step<(nsteps+1);++step){
    if(step==nsteps) fracstep=last_fracstep;

    unsigned last_position(0);
    unsigned this_position(0);
    unsigned next_position(1);

    // bootstrap the process by loading this_position before starting loop:
    n_ring[this_position] = static_cast<NUC::index_t>(leaf[0]);

    for(unsigned i(0);i<lsize;++i){
      // category at this position:
      const unsigned cat_index(cat_seq[i]);

      // push a new value in at next_position, if next_position exists:
      if((i+1)!=lsize){
        n_ring[next_position] = static_cast<NUC::index_t>(leaf[i+1]);
      }

      // For the first and last positions, choose the leading and
      // following nuc at random from the background distribution of
      // the first or last position category.
      //
      // This assumes that the random cdf draws, over many discrete
      // steps, have the same effect as averaging them:
      //
      NUC::index_t last_nuc(NUC::N);
      if       (i==0 || group_seq[i-1] != group_seq[i]){
        last_nuc = static_cast<NUC::index_t>(random_cdf_variate(bg_nuc_cdf_cat[cat_index],NUC::SIZE));
      } else { last_nuc=n_ring[last_position]; }

      NUC::index_t next_nuc(NUC::N);
      if(i==(lsize-1) || group_seq[i+1] != group_seq[i]) {
        next_nuc = static_cast<NUC::index_t>(random_cdf_variate(bg_nuc_cdf_cat[cat_index],NUC::SIZE));
      } else { next_nuc=n_ring[next_position]; }

      const unsigned nuc3mer(last_nuc+n_ring[this_position]*NUC::SIZE+next_nuc*(NUC::SIZE*NUC::SIZE));

      const double r(random_uniform());
      const prob_t* nucsub_cat(nucsub[cat_index]);
      if(r<(nucsub_cat[nuc3mer]*fracstep)){
        cat_sim_time[cat_index] += 1.;
        const prob_t* nucsub_cdf_cat(nucsub_cdf[cat_index]);
        int nuc = random_cdf_variate(nucsub_cdf_cat+nuc3mer*(NUC::SIZE-1),(NUC::SIZE-1));
        if(static_cast<int>(n_ring[this_position])<=nuc) nuc += 1;
        leaf[i] = nuc;
      }

      // increment positions
      last_position = this_position;
      this_position = next_position;
      next_position = (next_position+1)%RING_SIZE;
    }
  }
}



/// \todo add a codon-position specific bg_nuc distro
///
static
void
simulate_steps_codon(const std::vector<unsigned>& cat_seq,
                     const std::vector<unsigned>& group_seq,
                     std::vector<unsigned>& endseq,
                     const prob_t* const * codonsub,
                     const prob_t* const * codonsub_cdf,
                     const prob_t* const * bg_nuc_cdf_cat,
                     const bool is_report_sim_time,
                     const prob_t* const * codonsub_neutral,
                     const unsigned nsteps,
                     const smlfloat last_fracstep,
                     std::vector<smlfloat>& cat_sim_time){

  static const unsigned NUS(NUC::SIZE-1);

  // ring buffer (better described as a tiny fifo?) holds the +/-1 codon view
  // around the current codon
  static const unsigned RING_SIZE(3);

  NSCODON::index_t c_ring[RING_SIZE];
  NUC::index_t cn_ring[RING_SIZE][CODON::BASE_SIZE];

  const unsigned lsize(endseq.size());

  // precalc some indices:
  simple_array<sst_struct> pindex(lsize);
  for(unsigned i(0);i<lsize;++i){
    sst_struct& pi(pindex[i]);
    pi.cat=cat_seq[i];
    pi.break_5p=(i==0 || group_seq[i-1] != group_seq[i]);
    pi.break_3p=((i+1)==lsize || group_seq[i+1] != group_seq[i]);
  }

  const unsigned n_cats(cat_sim_time.size());
  for(unsigned i(0);i<n_cats;++i) cat_sim_time[i] = 0.;

  smlfloat fracstep(1.);
  NUC::index_t last_nuc(NUC::N);
  NUC::index_t next_nuc(NUC::N);

  for(unsigned step(0);step<(nsteps+1);++step){
    if(step==nsteps) fracstep=last_fracstep;

    unsigned last_position(0);
    unsigned this_position(0);
    unsigned next_position(1);

    // bootstrap the process by loading this_position before starting loop:
    c_ring[this_position] = static_cast<NSCODON::index_t>(endseq[0]);
    NSCODON::decode(cn_ring[this_position],c_ring[this_position]);

    for(unsigned i(0);i<lsize;++i){

      // push a new value in at next_position, unless at the last
      // position:
      if((i+1)!=lsize){
        c_ring[next_position] = static_cast<NSCODON::index_t>(endseq[i+1]);
        NSCODON::decode(cn_ring[next_position],c_ring[next_position]);
      }

      const sst_struct& pi(pindex[i]);
      const unsigned cat_index(pi.cat);

      // For the first and last positions, choose the leading and
      // following nuc at random from the background distribution of
      // the first or last position category.
      //
      // This assumes that the random cdf draws, over many discrete
      // steps, have the same effect as averaging them:
      //
      if(pi.break_5p){
        last_nuc = static_cast<NUC::index_t>(random_cdf_variate(bg_nuc_cdf_cat[cat_index],NUC::SIZE));
      } else { last_nuc=cn_ring[last_position][2]; }

      if(pi.break_3p) {
        next_nuc = static_cast<NUC::index_t>(random_cdf_variate(bg_nuc_cdf_cat[cat_index],NUC::SIZE));
      } else { next_nuc=cn_ring[next_position][0]; }

      const unsigned ns5mer = NSC4::encode(last_nuc,c_ring[this_position])+next_nuc*NSC4::SIZE;

      const double r(random_uniform());
      if(is_report_sim_time){
        const prob_t* codonsub_cat_neutral(codonsub_neutral[cat_index]);
        if(r<(codonsub_cat_neutral[ns5mer]*fracstep)) cat_sim_time[cat_index] += 1.;
      }

      const prob_t* codonsub_cat = codonsub[cat_index];
      if(r<(codonsub_cat[ns5mer]*fracstep)){
        const prob_t* codonsub_cdf_cat = codonsub_cdf[cat_index];
        const unsigned a = random_cdf_variate(codonsub_cdf_cat+ns5mer*CODON_CDF_SIZE,CODON_CDF_SIZE);
        const unsigned pos = a/NUS;
        unsigned nuc = a%NUS;
        if(static_cast<unsigned>(cn_ring[this_position][pos])<=nuc) nuc += 1;

        const NUC::index_t* nref = cn_ring[this_position];
        NUC::index_t ntmp[CODON::BASE_SIZE] = {nref[0],nref[1],nref[2]};
        ntmp[pos] = static_cast<NUC::index_t>(nuc);

        endseq[i] = NSCODON::encode(ntmp);
      }

      // increment positions
      last_position = this_position;
      this_position = next_position;
      next_position = (next_position+1)%RING_SIZE;
    }
  }
}



void
simulate_discrete_time_branch_context(const std::vector<unsigned>& startseq,
                                      const std::vector<unsigned>& cat_seq,
                                      const std::vector<unsigned>& group_seq,
                                      std::vector<unsigned>& endseq,
                                      const smlfloat time,
                                      const smlfloat unit_time,
                                      const prob_t* const * sitesub_prob,
                                      const prob_t* const * sitesub_cdf,
                                      const prob_t * const * bg_nuc_cdf_cat,
                                      std::vector<smlfloat>& cat_sim_time,
                                      const bool is_codon_site_model,
                                      const bool is_report_sim_time,
                                      const prob_t* const * sitesub_prob_neutral){

  const smlfloat tmp(time/unit_time);
  const unsigned nsteps(static_cast<unsigned>(tmp));
  const smlfloat fracstep(tmp-static_cast<smlfloat>(nsteps));

  log_os << "calc branch: steps,frac: " << nsteps << " " << fracstep << "\n";

  endseq = startseq;

  if(is_codon_site_model){
    simulate_steps_codon(cat_seq,group_seq,endseq,sitesub_prob,sitesub_cdf,
                         bg_nuc_cdf_cat,is_report_sim_time,sitesub_prob_neutral,
                         nsteps,fracstep,cat_sim_time);
  } else {
    simulate_steps_nuc(cat_seq,group_seq,endseq,sitesub_prob,sitesub_cdf,
                       bg_nuc_cdf_cat,nsteps,fracstep,cat_sim_time);
  }

  unsigned time_nuc_factor(1);
  if(is_codon_site_model) time_nuc_factor=CODON::BASE_SIZE;
  const smlfloat scale_factor(1./static_cast<smlfloat>(startseq.size()*time_nuc_factor));

  const unsigned n_cats(cat_sim_time.size());
  for(unsigned i(0); i<n_cats; ++i) cat_sim_time[i] *= scale_factor;
}

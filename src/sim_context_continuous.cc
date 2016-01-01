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

// $Id: sim_context_continuous.cc 1153 2008-03-18 00:24:07Z ctsa $

/// \file

#include "sim_context_continuous.h"
#include "sim_context_shared.h"
#include "simple_util.h"
#include "util/bio/bioseq_util.h"
#include "util/general/die.h"
#include "util/general/log.h"
#include "util/math/prob_util.h"
#include "util/math/random_util.h"

#include <cmath>

#include <ostream>



static
void
update_codon_pos_rates(const prob_t* bg_nuc_pdf,
                       const sst_struct& pi,
                       const prob_t* codonsub_cat,
                       const NSC5::index_t ns5mer,
                       smlfloat& codon_rate,
                       smlfloat& total_rate){

  NUC::index_t nsc5_nuc[SITE_MODEL::MAX_BASE_SIZE];

  if(pi.break_5p || pi.break_3p){
    // avg_pos = pos in ns5mer with 'unknown' nucleotide that we average over
    unsigned avg_pos=4;
    if(pi.break_5p) avg_pos=0;

    codon_rate=0;
    SITE_MODEL::decode_nuc(SITE_MODEL::NSC5,ns5mer,nsc5_nuc);
    for(unsigned n(0);n<NUC::SIZE;++n){
      nsc5_nuc[avg_pos]=static_cast<NUC::index_t>(n);
      const NSC5::index_t edge_ns5mer(SITE_MODEL::encode_nuc(SITE_MODEL::NSC5,nsc5_nuc));
      codon_rate += codonsub_cat[edge_ns5mer]*bg_nuc_pdf[n];
    }
  } else {
    codon_rate = codonsub_cat[ns5mer];
  }
  total_rate += codon_rate;
}



static
void
update_codon_pos(const prob_t* codonsub_cat,
                 const unsigned nsc5_pos,
                 const NUC::index_t nsc5_new_nuc,
                 const prob_t* bg_nuc_pdf,
                 const sst_struct& pi,
                 NSC5::index_t& codon_ns5mer,
                 smlfloat& codon_rate,
                 smlfloat& total_rate){

  NUC::index_t nsc5_nuc[SITE_MODEL::MAX_BASE_SIZE];

  SITE_MODEL::decode_nuc(SITE_MODEL::NSC5,codon_ns5mer,nsc5_nuc);
  nsc5_nuc[nsc5_pos]=nsc5_new_nuc;
  codon_ns5mer = SITE_MODEL::encode_nuc(SITE_MODEL::NSC5,nsc5_nuc);

  total_rate -= codon_rate;
  update_codon_pos_rates(bg_nuc_pdf,pi,codonsub_cat,codon_ns5mer,codon_rate,total_rate);
}


/// \todo the pointers here are only proxies for 'mutable', so break
/// this into two structs, one const and the other not
///
struct upcf_type {
  const prob_t* codonsub_cat;
  NUC::index_t nsc5_new_nuc;
  const prob_t* bg_nuc_pdf;
  const sst_struct* pi;
  NSC5::index_t* codon_ns5mer;
  smlfloat* codon_rate;
  smlfloat* total_rate_ptr;
  random_scaled_pdistro_variate_cached<smlfloat>* rspc_ptr;
  bool is_report_sim_time;
  const prob_t* codonsub_neutral_cat;
  smlfloat* codon_neutral_rate;
  std::vector<smlfloat>* total_neutral_rate_ptr;
};



static
void
update_codon_pos_full(const unsigned nsc5_pos,
                      const unsigned mut_codon,
                      const upcf_type& u){

  const smlfloat old_total_rate(*(u.total_rate_ptr));

  update_codon_pos(u.codonsub_cat,nsc5_pos,u.nsc5_new_nuc,
                   u.bg_nuc_pdf,u.pi[mut_codon],
                   u.codon_ns5mer[mut_codon],
                   u.codon_rate[mut_codon],
                   *(u.total_rate_ptr));

  u.rspc_ptr->update(mut_codon,*(u.total_rate_ptr)-old_total_rate);

  if(u.is_report_sim_time){
    smlfloat& neutral_rate(u.total_neutral_rate_ptr->operator[](u.pi[mut_codon].cat));
    neutral_rate -= u.codon_neutral_rate[mut_codon];
    update_codon_pos_rates(u.bg_nuc_pdf,u.pi[mut_codon],u.codonsub_neutral_cat,
                           u.codon_ns5mer[mut_codon],
                           u.codon_neutral_rate[mut_codon],
                           neutral_rate);
  }
}



/// \todo use codon position specific bg_nuc distro
///
static
void
simulate_branch_codon(const std::vector<unsigned>& cat_seq,
                      const std::vector<unsigned>& group_seq,
                      std::vector<unsigned>& endseq,
                      const prob_t* const * codonsub,
                      const prob_t* const * codonsub_cdf,
                      const prob_t* const * bg_nuc_cdf_cat,
                      const bool is_report_sim_time,
                      const prob_t* const * codonsub_neutral,
                      const smlfloat branch_time,
                      const unsigned n_cats,
                      std::vector<smlfloat>& cat_neutral_mutation){

  static const unsigned NUS(NUC::SIZE-1);

  // ring buffer (better described as a tiny fifo?) holds the +/-1 codon view
  // around the current codon
  static const unsigned RING_SIZE(3);

  NSCODON::index_t c_ring[RING_SIZE];
  NUC::index_t cn_ring[RING_SIZE][CODON::BASE_SIZE];

  const unsigned lsize(endseq.size());

  // precalc some indices that have no dependency on the state sequence:
  simple_array<sst_struct> pindex(lsize);
  for(unsigned i(0);i<lsize;++i){
    sst_struct& pi(pindex[i]);
    pi.cat=cat_seq[i];
    pi.break_5p=(i==0 || group_seq[i-1] != group_seq[i]);
    pi.break_3p=((i+1)==lsize || group_seq[i+1] != group_seq[i]);
  }

  NUC::index_t last_nuc(NUC::N);
  NUC::index_t next_nuc(NUC::N);

  // start:
  // get total sequence mutation rate
  // store total rate array
  //
  // per iteration:
  // ) find time to next mut
  // ) if out of time break
  //
  // ) pull random variate from sequence rate distribution to determine
  //   mutated codon
  // ) pull random variate from codon rate distribution to determine
  //   mutated nucleotide
  //
  // ) mutate state sequence
  // ) modify total sequence mutation rate
  // ) modify total rate array
  //

  // 1) get initial total sequence mutation rate & get sequence rate array:
  NUC::index_t nsc5_nuc[CODON::BASE_SIZE+2];

  smlfloat rate(0);
  simple_array<smlfloat> codon_rates(lsize);
  simple_array<NSC5::index_t> codon_ns5mer(lsize);

  std::vector<smlfloat> neutral_rate(n_cats,0);
  simple_init_array<smlfloat> codon_neutral_rates;
  if(is_report_sim_time){
    codon_neutral_rates.init(lsize);
  }

  simple_matrix<prob_t> bg_nuc_pdf_cat(n_cats,NUC::SIZE);
  for(unsigned c(0);c<n_cats;++c){
    pdistro_from_cdf(bg_nuc_cdf_cat[c],bg_nuc_pdf_cat[c],NUC::SIZE);
  }

  {
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

      if(pi.break_5p) { last_nuc = NUC::A; }
      else            { last_nuc=cn_ring[last_position][2]; }

      if(pi.break_3p) { next_nuc = NUC::A; }
      else            { next_nuc=cn_ring[next_position][0]; }

      const unsigned ns5mer(NSC4::encode(last_nuc,c_ring[this_position])+next_nuc*NSC4::SIZE);
      const unsigned& cat_index(pi.cat);

      codon_ns5mer[i] = ns5mer;
      update_codon_pos_rates(bg_nuc_pdf_cat[cat_index],pi,codonsub[cat_index],ns5mer,codon_rates[i],rate);
      if(is_report_sim_time){
        update_codon_pos_rates(bg_nuc_pdf_cat[cat_index],pi,codonsub_neutral[cat_index],ns5mer,
                               codon_neutral_rates[i],neutral_rate[cat_index]);
      }

      // increment positions
      last_position = this_position;
      this_position = next_position;
      next_position = (next_position+1)%RING_SIZE;
    }
  }

  assert(! is_float_invalid(rate));

  if(is_report_sim_time){
    for(unsigned i(0);i<n_cats;++i) cat_neutral_mutation[i] = 0.;
  }

  NUC::index_t codon_nuc[CODON::BASE_SIZE];

  bool is_break(false);
  smlfloat elapsed_time(0.);
  smlfloat n_check(0);
  unsigned n_check_i(0);
  random_scaled_pdistro_variate_cached<smlfloat> rspc(codon_rates.ptr(),lsize);
  while(true){
    // ) get time to next mutation:
    smlfloat mut_time(random_exponential(rate));

    if((elapsed_time+mut_time)>branch_time){
      mut_time=branch_time-elapsed_time;
      is_break=true;
    }

    elapsed_time += mut_time;

    if(is_report_sim_time) {
      for(unsigned c(0);c<n_cats;++c){
        cat_neutral_mutation[c] += neutral_rate[c]*mut_time;
      }
    }
    n_check += rate*mut_time;

    if(is_break) break;

    //tmp theoretical check
    //    log_os << "r,rate,mut_time: " << r << " " << rate << " " << mut_time << "\n";
    n_check_i++;

    // 3) pull cdf variate to determine which codon is mutated
    const unsigned mut_codon(rspc.get_random_draw());

    { // update nucleotide change to mutated codon:
      const sst_struct pi(pindex[mut_codon]);
      const unsigned cat_index(pi.cat);
      const NSC5::index_t ns5mer(codon_ns5mer[mut_codon]);
      const prob_t* codonsub_cdf_cat(codonsub_cdf[cat_index]);
      NSC5::index_t edge_ns5mer(ns5mer);
      if(pi.break_5p || pi.break_3p){
        unsigned edge_pos=4;
        if(pi.break_5p) edge_pos=0;

        const unsigned edge_nuc(random_cdf_variate(bg_nuc_cdf_cat[cat_index],NUC::SIZE));
        SITE_MODEL::decode_nuc(SITE_MODEL::NSC5,ns5mer,nsc5_nuc);
        nsc5_nuc[edge_pos]=static_cast<const NUC::index_t>(edge_nuc);
        edge_ns5mer=SITE_MODEL::encode_nuc(SITE_MODEL::NSC5,nsc5_nuc);
      }
      const unsigned codon_cdf_variate(random_cdf_variate(codonsub_cdf_cat+edge_ns5mer*CODON_CDF_SIZE,CODON_CDF_SIZE));
      const unsigned pos(codon_cdf_variate/NUS);
      unsigned nuc(codon_cdf_variate%NUS);

      NSCODON::decode(codon_nuc,static_cast<NSCODON::index_t>(endseq[mut_codon]));
      if(static_cast<unsigned>(codon_nuc[pos])<=nuc) nuc+=1;
      const NUC::index_t nuci(static_cast<NUC::index_t>(nuc));
      codon_nuc[pos]=nuci;

      endseq[mut_codon] = NSCODON::encode(codon_nuc);

      upcf_type ut;
      ut.codonsub_cat=codonsub[cat_index];
      ut.nsc5_new_nuc=nuci;
      ut.bg_nuc_pdf=bg_nuc_pdf_cat[cat_index];
      ut.pi=pindex.ptr();
      ut.codon_ns5mer=codon_ns5mer.ptr();
      ut.codon_rate=codon_rates.ptr();
      ut.total_rate_ptr=&rate;
      ut.rspc_ptr=&rspc;
      ut.is_report_sim_time=is_report_sim_time;
      if(is_report_sim_time){
        ut.codonsub_neutral_cat=codonsub_neutral[cat_index];
        ut.codon_neutral_rate=codon_neutral_rates.ptr();
        ut.total_neutral_rate_ptr=&neutral_rate;
      } else {
        ut.codonsub_neutral_cat=0;
        ut.codon_neutral_rate=0;
        ut.total_neutral_rate_ptr=0;
      }
      update_codon_pos_full(pos+1,mut_codon,ut);

      if       (pos==0 && mut_codon!=0){          //update previous codon:
        update_codon_pos_full(4,mut_codon-1,ut);
      } else if(pos==2 && (mut_codon+1)!=lsize){  //update next codon:
        update_codon_pos_full(0,mut_codon+1,ut);
      }
    }
  }

  log_os << "nchk,nchki: " << n_check << " " << n_check_i << " " << n_check/n_check_i << "\n";
}




void
simulate_continuous_time_branch_context(const std::vector<unsigned>& startseq,
                                        const std::vector<unsigned>& cat_seq,
                                        const std::vector<unsigned>& group_seq,
                                        std::vector<unsigned>& endseq,
                                        const smlfloat time,
                                        const unsigned n_cats,
                                        const prob_t* const * sitesub_prob,
                                        const prob_t* const * sitesub_cdf,
                                        const prob_t* const * bg_nuc_cdf_cat,
                                        std::vector<smlfloat>& cat_sim_time,
                                        const bool is_codon_site_model,
                                        const bool is_report_sim_time,
                                        const prob_t* const * sitesub_prob_neutral){

  log_os << "calc branch: time: " << time << "\n";

  endseq = startseq;

  if(is_codon_site_model){
    simulate_branch_codon(cat_seq,group_seq,endseq,sitesub_prob,sitesub_cdf,
                          bg_nuc_cdf_cat,is_report_sim_time,sitesub_prob_neutral,
                          time,n_cats,cat_sim_time);
  } else {
    pass_away("noncoding model not yet supported in continuous time simulator");
    //    n_mutation = simulate_steps_nuc(site_cats,group_cats,group_seq,endseq,sitesub_prob,sitesub_cdf,
    //                                    bg_nuc_cdf_cat,n_site_cats,nsteps,fracstep);
  }

  unsigned time_nuc_factor(1);
  if(is_codon_site_model) time_nuc_factor=CODON::BASE_SIZE;

  if(is_report_sim_time){
    const smlfloat time_scale(1./static_cast<smlfloat>(startseq.size()*time_nuc_factor));

    for(unsigned i(0);i<n_cats;++i) cat_sim_time[i] *= time_scale;
  }
}

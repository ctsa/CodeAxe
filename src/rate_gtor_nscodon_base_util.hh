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

// $Id: rate_gtor_nscodon_base_util.hh 1153 2008-03-18 00:24:07Z ctsa $

/// \file

#include "context_group.h"
#include "simple_util.h"
#include "util/bio/bioseq_util_pdf_conversions.h"
#include "util/bio/bioseq_util_flux_rate_conversions.h"
#include "util/math/array_util.h"


template <typename FloatType>
void
nscodon_flux_no_aa_selection(const rate_gtor_nscodon_base& r,
                             const unsigned cat_no,
                             const unsigned branch_cat_set,
                             FloatType* nscodon_flux) {

  const unsigned n_states(r.state_size());
  const unsigned N2(n_states*n_states);
  simple_array<FloatType> flux(N2,0.);

  r.rates_variance(flux.ptr(),rates_func_options(rates_func_options_base(cat_no,true),branch_cat_set));
  rates_to_flux_inplace(n_states,r.bg_pdistro_cat(cat_no),flux.ptr());
  site_model_flux_to_nscodon_flux(r.site_model(),flux.ptr(),nscodon_flux);
}



template <typename FloatType>
void
rates_nscodon_context(const rate_gtor_nscodon_base& r,
                      const SITE_MODEL::index_t sm,
                      FloatType rates[],
                      const rates_func_options& opt){

  const unsigned base_size(SITE_MODEL::base_size(sm));
  unsigned codon_pos[SITE_MODEL::MAX_BASE_SIZE];
  SITE_MODEL::codon_position(sm,codon_pos);

  const unsigned overedge_pos_pre(periodic_decrement(codon_pos[0],CODON::BASE_SIZE));
  const unsigned overedge_pos_post(periodic_increment(codon_pos[base_size-1],CODON::BASE_SIZE));

  const prob_t* nuc_pdistro_overedge_5p(opt.nuc_pdistro_overedge_5p ? opt.nuc_pdistro_overedge_5p :
                                        r.bg_pdistro_nuc_pos_cat(opt.cat,overedge_pos_pre));
  const prob_t* nuc_pdistro_overedge_3p(opt.nuc_pdistro_overedge_3p ? opt.nuc_pdistro_overedge_3p :
                                        r.bg_pdistro_nuc_pos_cat(opt.cat,overedge_pos_post));

  const prob_t* nuc_pdistro_overedge_5p_cond_on_edge(opt.nuc_pdistro_overedge_5p_cond_on_edge ? opt.nuc_pdistro_overedge_5p_cond_on_edge :
                                                     r.bg_pdistro_nuc_pos_conditioned_on_3p_cat(opt.cat,overedge_pos_pre));
  const prob_t* nuc_pdistro_overedge_3p_cond_on_edge(opt.nuc_pdistro_overedge_3p_cond_on_edge ? opt.nuc_pdistro_overedge_3p_cond_on_edge :
                                                     r.bg_pdistro_nuc_pos_conditioned_on_5p_cat(opt.cat,overedge_pos_post));

  const prob_t* nscodon_pdistro(opt.nscodon_pdistro ? opt.nscodon_pdistro :
                                r.bg_pdistro_nscodon_cat(opt.cat));


  prob_t overedge_pos_pre_cond_3p[NUC::SIZE][NUC::SIZE];
  prob_t overedge_pos_post_cond_5p[NUC::SIZE][NUC::SIZE];
  for(unsigned i(0);i<NUC::SIZE;++i){
    const NUC::index_t ni(static_cast<NUC::index_t>(i));
    for(unsigned j(0);j<NUC::SIZE;++j){
      overedge_pos_pre_cond_3p[i][j] =
        nuc_pdistro_overedge_5p_cond_on_edge[ni*NUC::SIZE+j]*r.edge_strength_pre()+
        nuc_pdistro_overedge_5p[j]*(1.-r.edge_strength_pre());

      overedge_pos_post_cond_5p[i][j] =
        nuc_pdistro_overedge_3p_cond_on_edge[ni*NUC::SIZE+j]*r.edge_strength_post()+
        nuc_pdistro_overedge_3p[j]*(1.-r.edge_strength_post());
    }
  }

  const unsigned state_size(SITE_MODEL::state_size(sm));

  bool is_partial_codon[SITE_MODEL::MAX_BASE_SIZE];
  SITE_MODEL::is_partial_codon(sm,is_partial_codon);

  NUC::index_t s1n[SITE_MODEL::MAX_BASE_SIZE];
  NUC::index_t s2n[SITE_MODEL::MAX_BASE_SIZE];

  // calculate actual rate matrix
  for(unsigned s1(0);s1<state_size;++s1){
    SITE_MODEL::decode_nuc(sm,s1,s1n);
    const NSCODON::index_t s1nc(SITE_MODEL::decode_nscodon(sm,s1));
    const NSAA::index_t s1aa(codon_trans_known(s1nc));

    for(unsigned s2(0);s2<state_size;++s2){
      if( s1==s2 ) continue;

      SITE_MODEL::decode_nuc(sm,s2,s2n);
      const NSCODON::index_t s2nc(SITE_MODEL::decode_nscodon(sm,s2));
      const NSAA::index_t s2aa(codon_trans_known(s2nc));

      unsigned nsubs(0);
      unsigned mut_p(0);

      for(unsigned p(0);p<base_size;++p){
        if( s1n[p] != s2n[p] ) {
          mut_p=p;
          nsubs++;
        }
      }

      irv_t<smlfloat> this_rate(0.);

      if( nsubs == 1 ) {
        const unsigned mut_codon_pos(codon_pos[mut_p]);
        const NUC::index_t s1ndiff(s1n[mut_p]);
        const NUC::index_t s2ndiff(s2n[mut_p]);
        const NUC::index_t s1ndiff_plus1 ( ((mut_p+1)<base_size) ? s1n[mut_p+1] : NUC::N );
        const NUC::index_t s1ndiff_minus1( (mut_p>0)             ? s1n[mut_p-1] : NUC::N );

        const prob_t* bg_minus1(0);
        if(s1ndiff_minus1 == NUC::N) { bg_minus1 = overedge_pos_pre_cond_3p[s1ndiff]; }

        const prob_t* bg_plus1(0);
        if(s1ndiff_plus1 == NUC::N) { bg_plus1 = overedge_pos_post_cond_5p[s1ndiff]; }

        this_rate = r.cat_nuc_context_mut_rate(s1ndiff_minus1,s1ndiff,s1ndiff_plus1,s2ndiff,
                                               opt.cat,opt.branch_cat_set,bg_minus1,bg_plus1,s2nc);

        if( ! is_partial_codon[mut_p]){
          // selection when full codon context is known:
          //

          // codon bias selection
          if(! opt.is_skip_nonaa_selection){
            this_rate *= r.codon_bias(s1nc,s2nc);
          }

          // non-synonymous aa exchange selection:
          if( s1aa != s2aa && (! opt.is_skip_aa_selection)){
            this_rate *= r.cat_nsaa_param(s1aa,s2aa,opt.cat,opt.branch_cat_set);
          }
        } else {
          // selection when partial codon context is known:
          //
          // note that the skip_*aa_selection flags are accounted for
          // within average_codon_selection(); even for the all-skip
          // case, we still need to account for selection against stop
          // codons:
          //

          /// \todo get this to work for multiple (2) nucs fixed in partial codon
          NUC::index_t nuc[NUC::SIZE];
          for(unsigned i(0);i<CODON::BASE_SIZE;++i){
            if(i==mut_codon_pos) { nuc[i] = s1ndiff; }
            else                 { nuc[i] = NUC::N; }
          }

          const smlfloat s(average_codon_selection(r,nuc,s2ndiff,mut_codon_pos,opt,
                                                   nscodon_pdistro));
          this_rate *= s;

        }
      }

      rates[s2+s1*state_size] = this_rate;
    }
  }

  matrix_balance_diagonal(rates,state_size);
}

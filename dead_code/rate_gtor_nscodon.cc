// -*- mode: c++; indent-tabs-mode: nil; -*-
//

/// \file rate_gtor_nscodon.cc
/// \brief implementation for non-stop codon model rate generator
///


#include "context_group.h"
#include "rate_edge_dependencies.h"
#include "rate_gtor_nscodon.h"
#include "util/math/matrix_util.h"
#include "rate_gtor_nscodon_base_util.h"

#include <cassert>

#include <iostream>


// put rate information into a c-style matrix
void
rate_gtor_nscodon::
rates(smlfloat rates[],
      const rates_func_options& opt) const {

  rates_nscodon_context(*this,rates,opt);
#if 0
  const CONTEXT_GROUP::index_t cg(get_context_model_group(context_model_nuc()));


  // first calculate edge dependencies:
  prob_t bg_pdf_nuc0_cond_on_nuc2[NUC::SIZE*NUC::SIZE];
  prob_t bg_pdf_nuc2_cond_on_nuc0[NUC::SIZE*NUC::SIZE];
  if(opt.is_use_edge_correction){
    if(cg == CONTEXT_GROUP::ONLY_POST || cg == CONTEXT_GROUP::TWOWAY){
      static const int edge_pos(0);
      static const int overedge_dir(-1);
      get_rate_edge_dependencies_coding(*this,opt,edge_pos,overedge_dir,
                                        0,bg_pdf_nuc0_cond_on_nuc2);
    }
    if(cg == CONTEXT_GROUP::ONLY_PRE || cg == CONTEXT_GROUP::TWOWAY){
      static const int edge_pos(2);
      static const int overedge_dir(1);
      get_rate_edge_dependencies_coding(*this,opt,edge_pos,overedge_dir,
                                        0,bg_pdf_nuc2_cond_on_nuc0);
    }
  }

  // now build the primary rate matrix:
  for(unsigned c1(0);c1<NSCODON::SIZE;++c1){
    NUC::index_t c1n1,c1n2,c1n3;
    NSCODON::decode(c1n1,c1n2,c1n3,static_cast<NSCODON::index_t>(c1));
    NSAA::index_t c1aa = codon_trans_known(static_cast<NSCODON::index_t>(c1));

    for(unsigned c2(0);c2<NSCODON::SIZE;++c2){
      if( c1==c2 ) continue;
      NUC::index_t c2n1,c2n2,c2n3;
      NSCODON::decode(c2n1,c2n2,c2n3,static_cast<NSCODON::index_t>(c2));
      NSAA::index_t c2aa = codon_trans_known(static_cast<NSCODON::index_t>(c2));

      unsigned nsubs(0);
      NUC::index_t c1ndiff(NUC::N);
      NUC::index_t c1ndiff_minus1(NUC::N);
      NUC::index_t c1ndiff_plus1(NUC::N);
      NUC::index_t c2ndiff(NUC::N);
      unsigned pos_mut(0);

      if( c1n1 != c2n1 ) {
        pos_mut = 0;
        nsubs++; c1ndiff=c1n1; c2ndiff=c2n1;
        c1ndiff_plus1  = c1n2;
      }
      if( c1n2 != c2n2 ) {
        pos_mut = 1;
        nsubs++; c1ndiff=c1n2; c2ndiff=c2n2;
        c1ndiff_minus1 = c1n1;
        c1ndiff_plus1  = c1n3;
      }
      if( c1n3 != c2n3 ) {
        pos_mut = 2;
        nsubs++; c1ndiff=c1n3; c2ndiff=c2n3;
        c1ndiff_minus1 = c1n2;
      }

      smlfloat this_rate(0.);

      if( nsubs == 1 ) {
        const int pos_minus1 = (pos_mut+(CODON::BASE_SIZE-1))%CODON::BASE_SIZE;
        const int pos_plus1 = (pos_mut+1)%CODON::BASE_SIZE;

        const prob_t* bg_minus1 = bg_pdistro_nuc_pos_opt(pos_minus1);
        const prob_t* bg_plus1 = bg_pdistro_nuc_pos_opt(pos_plus1);
        if( opt.is_use_edge_correction ){
          if((cg == CONTEXT_GROUP::ONLY_PRE || cg == CONTEXT_GROUP::TWOWAY) && pos_minus1==2){
            bg_minus1 = bg_pdf_nuc2_cond_on_nuc0+c1ndiff*NUC::SIZE;
          }
          if((cg == CONTEXT_GROUP::ONLY_POST || cg == CONTEXT_GROUP::TWOWAY) && pos_plus1==0){
            bg_plus1 = bg_pdf_nuc0_cond_on_nuc2+c1ndiff*NUC::SIZE;
          }
        }

        // mutation rate
        this_rate = get_category_nuc_param(c1ndiff_minus1,c1ndiff,c1ndiff_plus1,c2ndiff,
                                           opt.site_category,opt.group_category,
                                           bg_minus1,bg_plus1);

        // codon bias..
        if( ! opt.is_skip_nonaa_selection ){
          this_rate *= get_codon_bias(static_cast<NSCODON::index_t>(c1),static_cast<NSCODON::index_t>(c2));
        }
        if( ( c1aa != c2aa ) && ( ! opt.is_skip_aa_selection ) ) {
          this_rate *= get_category_nsaa_param(c1aa,c2aa,opt.site_category,opt.group_category);
        }
      }

      // finished.. store rates
      rates[c2+c1*NSCODON::SIZE] = this_rate;
    }
  }

  matrix_balance_diagonal(rates,NSCODON::SIZE);
#endif
}

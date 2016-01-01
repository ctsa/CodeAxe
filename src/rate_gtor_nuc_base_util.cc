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

// $Id: rate_gtor_nuc_base_util.cc 966 2007-10-26 17:42:16Z ctsa $

/// \file

//#include "rate_edge_dependencies.h"
#include "rate_gtor_nuc_base_util.h"
#include "util/bio/bioseq_util_pdf_conversions.h"
#include "util/math/matrix_util.h"


// nuc model rate matrix \w context
void
rates_nuc_context(const rate_gtor_nuc_base& r,
                  smlfloat rates[],
                  const rates_func_options& opt) {

  const prob_t* nuc_pdistro_overedge_5p(opt.nuc_pdistro_overedge_5p ? opt.nuc_pdistro_overedge_5p :
                                        r.bg_pdistro_nuc_cat(opt.cat));
  const prob_t* nuc_pdistro_overedge_3p(opt.nuc_pdistro_overedge_3p ? opt.nuc_pdistro_overedge_3p :
                                        r.bg_pdistro_nuc_cat(opt.cat));
  const prob_t* nuc_pdistro_overedge_5p_cond_on_edge(opt.nuc_pdistro_overedge_5p_cond_on_edge ? opt.nuc_pdistro_overedge_5p_cond_on_edge :
                                                     r.bg_pdistro_nuc_conditioned_on_3p_cat(opt.cat));
  const prob_t* nuc_pdistro_overedge_3p_cond_on_edge(opt.nuc_pdistro_overedge_3p_cond_on_edge ? opt.nuc_pdistro_overedge_3p_cond_on_edge :
                                                     r.bg_pdistro_nuc_conditioned_on_5p_cat(opt.cat));

  prob_t overedge_pre_cond_3p[NUC::SIZE][NUC::SIZE];
  prob_t overedge_post_cond_5p[NUC::SIZE][NUC::SIZE];
  for(unsigned i(0);i<NUC::SIZE;++i){
    const NUC::index_t ni(static_cast<NUC::index_t>(i));
    for(unsigned j(0);j<NUC::SIZE;++j){
      overedge_pre_cond_3p[i][j] =
        nuc_pdistro_overedge_5p_cond_on_edge[ni*NUC::SIZE+j]*r.edge_strength_pre()+
        nuc_pdistro_overedge_5p[j]*(1.-r.edge_strength_pre());

      overedge_post_cond_5p[i][j] =
        nuc_pdistro_overedge_3p_cond_on_edge[ni*NUC::SIZE+j]*r.edge_strength_post()+
        nuc_pdistro_overedge_3p[j]*(1.-r.edge_strength_post());
    }
  }


  const SITE_MODEL::index_t sm(r.site_model());

  const unsigned base_size(SITE_MODEL::base_size(sm));
  const unsigned state_size(SITE_MODEL::state_size(sm));

  NUC::index_t s1n[SITE_MODEL::MAX_BASE_SIZE];
  NUC::index_t s2n[SITE_MODEL::MAX_BASE_SIZE];

  for(unsigned s1(0);s1<state_size;++s1){
    SITE_MODEL::decode_nuc(sm,s1,s1n);

    for(unsigned s2(0);s2<state_size;++s2){
      if( s1==s2 ) continue;
      SITE_MODEL::decode_nuc(sm,s2,s2n);

      unsigned nsubs(0);
      unsigned mut_p(0);

      for(unsigned p(0);p<base_size;++p){
        if( s1n[p] != s2n[p] ){
          nsubs++;
          mut_p=p;
        }
      }

      smlfloat this_rate(0.);

      if( nsubs == 1 ) {
        const NUC::index_t s1ndiff(s1n[mut_p]);
        const NUC::index_t s2ndiff(s2n[mut_p]);
        const NUC::index_t s1ndiff_plus1 ( ((mut_p+1)<base_size) ? s1n[mut_p+1] : NUC::N );
        const NUC::index_t s1ndiff_minus1( (mut_p>0)             ? s1n[mut_p-1] : NUC::N );

        const prob_t* bg_minus1(0);
        if(s1ndiff_minus1 == NUC::N) { bg_minus1 = overedge_pre_cond_3p[s1ndiff]; }

        const prob_t* bg_plus1(0);
        if(s1ndiff_plus1 == NUC::N) { bg_plus1 = overedge_post_cond_5p[s1ndiff]; }


        this_rate = r.cat_nuc_context_mut_rate(s1ndiff_minus1,s1ndiff,s1ndiff_plus1,s2ndiff,
                                               opt.cat,opt.branch_cat_set,bg_minus1,bg_plus1);
      }
      // finished.. store rates
      rates[s2+s1*state_size] = this_rate;
    }
  }

  matrix_balance_diagonal(rates,state_size);
}

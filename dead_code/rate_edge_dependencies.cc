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

// $Id: rate_edge_dependencies.cc 1099 2008-01-21 22:54:47Z ctsa $

/// \file

#include "rate_edge_dependencies.h"
#include "simple_util.h"
#include "util/bio/bioseq_util_pdf_conversions.h"
#include "util/general/log.h"


#include <cmath>
#include <ostream>

#if 0

// manipulate joint prob to have assigned marginal probabilities
//
static
void
shape_joint_nuc_pdistro_margins(const prob_t* bg_pdf_pos1,
                                const prob_t* bg_pdf_pos2,
                                prob_t p_joint[NUC::SIZE][NUC::SIZE]){

  static const int MAXITER(100);
  for(int iter(0);true;iter++){
    if(iter>MAXITER){ die("can't shape joint nuc distro\n"); }

    // p(1,2) to p(1|2):
    for(int i(0);i<NUC::SIZE;++i){
      prob_t sum(0.);
      for(int j(0);j<NUC::SIZE;++j){ sum += p_joint[j][i]; }
      for(int j(0);j<NUC::SIZE;++j){ p_joint[j][i] /= sum; }
    }

    // test p(1)
    bool p_1_converge(true);
    for(int j(0);j<NUC::SIZE;++j){
      prob_t sum(0.);
      for(int i(0);i<NUC::SIZE;++i){
        sum += p_joint[j][i]*bg_pdf_pos2[i];
      }
      if(std::fabs(sum-bg_pdf_pos1[j])>0.001){
        p_1_converge=false;
        break;
      }
    }

    // p(1|2) to p(1,2)
    for(int i(0);i<NUC::SIZE;++i){
      for(int j(0);j<NUC::SIZE;++j){
        p_joint[i][j] *= bg_pdf_pos2[j];
      }
    }

    // p(1,2) to p(2|1)
    for(int i(0);i<NUC::SIZE;++i){
      prob_t sum(0.);
      for(int j(0);j<NUC::SIZE;++j){ sum += p_joint[i][j]; }
      for(int j(0);j<NUC::SIZE;++j){ p_joint[i][j] /= sum; }
    }

    // test p(2)
    bool p_2_converge(true);
    for(int i(0);i<NUC::SIZE;++i){
      prob_t sum(0.);
      for(int j(0);j<NUC::SIZE;++j){
        sum += p_joint[j][i]*bg_pdf_pos1[j];
      }
      if(std::fabs(sum-bg_pdf_pos2[i])>0.001){
        p_2_converge=false;
        break;
      }
    }

    // p(2|1) to p(1,2)
    for(int i(0);i<NUC::SIZE;++i){
      for(int j(0);j<NUC::SIZE;++j){
        p_joint[i][j] *= bg_pdf_pos1[i];
      }
    }

    if(p_1_converge && p_2_converge) return;
  }
}


// manipulate dinuc distribution to have marginal
// probabilities of pos 2,0 from nscodon distro
//
void
shape_dinuc_nscodon(const prob_t* x_nscodon,
                    prob_t* x_dinuc){

  prob_t x_dinuc2[NUC::SIZE][NUC::SIZE];
  for(unsigned i(0);i<DINUC::SIZE;++i){
    const NUC::index_t nx = DINUC::decode_nx(i);
    const NUC::index_t n0 = DINUC::decode_n0(i);
    x_dinuc2[nx][n0] = x_dinuc[i];
  }
  {
    prob_t n0[NUC::SIZE];
    nscodon_pdf_2_nuc_pdf_pos(x_nscodon,0,n0);
    prob_t n2[NUC::SIZE];
    nscodon_pdf_2_nuc_pdf_pos(x_nscodon,2,n2);

    shape_joint_nuc_pdistro_margins(n2,n0,x_dinuc2);
  }
  for(unsigned i(0);i<DINUC::SIZE;++i){
    const NUC::index_t nx = DINUC::decode_nx(i);
    const NUC::index_t n0 = DINUC::decode_n0(i);
    x_dinuc[i] = x_dinuc2[nx][n0];
  }
}

#endif


#ifdef BRANCH_SPECIFIC_OBS
void
site_model_edge_pdistros::
ropt_hookup(rates_func_options& ropt) const {

  ropt.nuc_pdistro_overedge_5p = nuc_pdistro_overedge_5p;
  ropt.nuc_pdistro_overedge_3p = nuc_pdistro_overedge_3p;
  ropt.nuc_pdistro_overedge_5p_cond_on_edge = nuc_pdistro_overedge_5p_cond_on_edge;
  ropt.nuc_pdistro_overedge_3p_cond_on_edge = nuc_pdistro_overedge_3p_cond_on_edge;
  ropt.nscodon_pdistro = nscodon_pdistro;
}



void
set_site_model_edge_pdistros(const SITE_MODEL::index_t sm,
                             const prob_t* sm_pdistro,
                             site_model_edge_pdistros& smp){
  if(sm != SITE_MODEL::NSC4PRE && sm != SITE_MODEL::NSC4POST){
    die("experimental edge distro function not yet generalized to all site models.");
  }
  const unsigned base_size(SITE_MODEL::base_size(sm));
  unsigned codon_pos[SITE_MODEL::MAX_BASE_SIZE];
  SITE_MODEL::codon_position(sm,codon_pos);

  const unsigned overedge_pos_pre(periodic_decrement(codon_pos[0],CODON::BASE_SIZE));
  const unsigned overedge_pos_post(periodic_increment(codon_pos[base_size-1],CODON::BASE_SIZE));

  unsigned base_no(0);
  {
    bool is_base_no_3p(false);
    unsigned base_no_3p(0);

    for(unsigned b(0);b<base_size;++b){
      if(codon_pos[b] == overedge_pos_pre){
        base_no=b;
        if((b+1)<base_size){
          base_no_3p=b+1;
          is_base_no_3p=true;
          break;
        }
      }
    }
    get_nuc_pos_distro_from_site_model_distro(sm,sm_pdistro,base_no,smp.nuc_pdistro_overedge_5p);
    get_dependent_nuc_pos_distro_from_site_model_distro(sm,sm_pdistro,base_no,is_base_no_3p,base_no_3p,
                                                        smp. nuc_pdistro_overedge_5p_cond_on_edge);
  }

  {
    bool is_base_no_5p(false);
    unsigned base_no_5p(0);

    for(unsigned b(0);b<base_size;++b){
      if(codon_pos[b] == overedge_pos_post){
        base_no=b;
        if(b>0){
          base_no_5p=b-1;
          is_base_no_5p=true;
          break;
        }
      }
    }

    get_nuc_pos_distro_from_site_model_distro(sm,sm_pdistro,base_no,smp.nuc_pdistro_overedge_3p);
    get_dependent_nuc_pos_distro_from_site_model_distro(sm,sm_pdistro,base_no,is_base_no_5p,base_no_5p,
                                                        smp. nuc_pdistro_overedge_3p_cond_on_edge);
  }

  // get nscodon:
  const BIO_PDISTRO::index_t i(BIO_PDISTRO::convert_from_site_model(sm));
  bio_pdistro_convert(i,BIO_PDISTRO::NSCODON,sm_pdistro,smp.nscodon_pdistro);
}
#endif




#if 0

#include "stationary_pdistro.h"
#include "rate_gtor_nscodon_base_util.h"
#include "subs_ml_ptol.h"
#include "util/general/die.h"
#include "util/math/matrix_util.h"

#include <cstdlib>

#include <iostream>


static
void
finish_edge_dependency_distro(const prob_t* bg_pdf_edge_pos,
                              const prob_t* bg_pdf_overedge_pos,
                              prob_t p_edge_overedge_tmp[NUC::SIZE][NUC::SIZE],
                              prob_t* bg_pdf_nuc_overedge_cond_on_edge,
                              prob_t* bg_pdf_nuc_edge_cond_on_overedge){

  // p(oe|e) to p(oe,e)
  for(int i(0);i<NUC::SIZE;++i){
    for(int j(0);j<NUC::SIZE;++j){
      p_edge_overedge_tmp[i][j] *= bg_pdf_edge_pos[i];
    }
  }

#ifndef SKIP_SHAPING
  shape_joint_nuc_pdistro_margins(bg_pdf_edge_pos,bg_pdf_overedge_pos,p_edge_overedge_tmp);
#endif

  // p(oe,e) to p(oe|e)
  for(int i(0);i<NUC::SIZE;++i){
    prob_t sum(0.);
    for(int j(0);j<NUC::SIZE;++j){ sum += p_edge_overedge_tmp[i][j]; }
    for(int j(0);j<NUC::SIZE;++j){ p_edge_overedge_tmp[i][j] /= sum; }
  }

  if(bg_pdf_nuc_overedge_cond_on_edge){
#if 0
    log_os << "p(oe):\n";
    for(int j(0);j<NUC::SIZE;++j) log_os << NUC::syms[j] << " " << bg_pdf_overedge_pos[j] << "\n";
    log_os << "p(oe|e)\n";
    for(int i(0);i<NUC::SIZE;++i)
      for(int j(0);j<NUC::SIZE;++j) log_os << NUC::syms[j] << " " << NUC::syms[i] << " " << p_edge_overedge_tmp[j][i] << "\n";
#endif
    // store p(oe|e)
    for(int i(0);i<NUC::SIZE;++i){
      for(int j(0);j<NUC::SIZE;++j){
        bg_pdf_nuc_overedge_cond_on_edge[j+i*NUC::SIZE] = p_edge_overedge_tmp[i][j];
      }
    }
#ifdef DEBUG
    for(int i(0);i<NUC::SIZE;++i){
      pdistro_check(bg_pdf_nuc_overedge_cond_on_edge+i*NUC::SIZE,NUC::SIZE,SUBS_ML_PTOL);
    }
#endif
  }


  if(bg_pdf_nuc_edge_cond_on_overedge){
    // p(oe|e) to p(oe,e)
    for(int i(0);i<NUC::SIZE;++i){
      for(int j(0);j<NUC::SIZE;++j){
        p_edge_overedge_tmp[i][j] *= bg_pdf_edge_pos[i];
      }
    }

    // p(oe,e) to p(e|oe):
    for(int i(0);i<NUC::SIZE;++i){
      prob_t sum(0.);
      for(int j(0);j<NUC::SIZE;++j){ sum += p_edge_overedge_tmp[j][i]; }
      for(int j(0);j<NUC::SIZE;++j){ p_edge_overedge_tmp[j][i] /= sum; }
    }
#if 0
    log_os << "p(e):\n";
    for(int j(0);j<NUC::SIZE;++j) log_os << NUC::syms[j] << " " << bg_pdf_edge_pos[j] << "\n";
    log_os << "p(e|oe)\n";
    for(int i(0);i<NUC::SIZE;++i)
      for(int j(0);j<NUC::SIZE;++j) log_os << NUC::syms[i] << " " << NUC::syms[j] << " " << p_edge_overedge_tmp[i][j] << "\n";
#endif
    // store p(e|oe):
    for(int i(0);i<NUC::SIZE;++i){
      for(int j(0);j<NUC::SIZE;++j){
        bg_pdf_nuc_edge_cond_on_overedge[i+j*NUC::SIZE] = p_edge_overedge_tmp[i][j];
      }
    }
#ifdef DEBUG
    for(int i(0);i<NUC::SIZE;++i){
      pdistro_check(bg_pdf_nuc_edge_cond_on_overedge+i*NUC::SIZE,NUC::SIZE,SUBS_ML_PTOL);
    }
#endif
  }
}




static
void
finish_edge_dependency_distro_bidirectional_context(const prob_t* bg_pdf_edge_pos,
                                                    const prob_t* bg_pdf_overedge_pos,
                                                    prob_t p_edge_tmp[NUC::SIZE][NUC::SIZE][NUC::SIZE],
                                                    prob_t bg_pdf_nuc_overedge_cond_on_edge_underedge[NUC::SIZE][NUC::SIZE][NUC::SIZE]){

  typedef prob_t p_subedge_t[NUC::SIZE][NUC::SIZE];

  for(unsigned n_underedge(0);n_underedge<NUC::SIZE;++n_underedge){
    p_subedge_t& p_edge_overedge_tmp(p_edge_tmp[n_underedge]);

#ifndef SKIP_SHAPING
    shape_joint_nuc_pdistro_margins(bg_pdf_edge_pos,bg_pdf_overedge_pos,p_edge_overedge_tmp);
#endif

    // p(e,oe) to p(oe|e):
    for(int i(0);i<NUC::SIZE;++i){
      prob_t sum(0.);
      for(int j(0);j<NUC::SIZE;++j){ sum += p_edge_overedge_tmp[i][j]; }
      for(int j(0);j<NUC::SIZE;++j){ p_edge_overedge_tmp[i][j] /= sum; }
    }
#if 0
    log_os << "u= " << NUC::syms[n_underedge] << " p(e):\n";
    for(int j(0);j<NUC::SIZE;++j) log_os << NUC::syms[j] << " " << bg_pdf_edge_pos[j] << "\n";
    log_os << "u= " << NUC::syms[n_underedge] << " p(e|oe):\n";
    for(int i(0);i<NUC::SIZE;++i)
      for(int j(0);j<NUC::SIZE;++j) log_os << NUC::syms[i] << " " << NUC::syms[j] << " " << p_edge_overedge_tmp[i][j] << "\n";
#endif
    // store p(oe|ue,e):
    for(int i(0);i<NUC::SIZE;++i){
      //      DINUC::index_t d_ue_e = DINUC::encode(static_cast<NUC::index_t>(n_underedge),static_cast<NUC::index_t>(i));
      for(int j(0);j<NUC::SIZE;++j){
        bg_pdf_nuc_overedge_cond_on_edge_underedge[n_underedge][i][j] = p_edge_overedge_tmp[i][j];
      }
    }
  }
}




void
get_rate_edge_dependencies_noncoding_unidirectional_context(const rate_gtor_nuc_base& r,
                                                            const unsigned site_category,
                                                            const unsigned group_category,
                                                            const int overedge_dir,
                                                            prob_t* bg_pdf_nuc_overedge_cond_on_edge,
                                                            prob_t* bg_pdf_nuc_edge_cond_on_overedge){

  const bool is_edge_before_overedge(overedge_dir>0);

  smlfloat edge_tmp;
  if(is_edge_before_overedge){ edge_tmp=r.edge_strength_pre(); }
  else                       { edge_tmp=r.edge_strength_post(); }

  const smlfloat edge_strength(edge_tmp);

  prob_t p_edge_overedge_tmp[NUC::SIZE][NUC::SIZE];

  for(unsigned n_edge(0);n_edge<NUC::SIZE;++n_edge){

    const prob_t* bg_pdf_pos_before_overedge = r.bg_pdistro_nuc();
    const prob_t* bg_pdf_pos_after_overedge = r.bg_pdistro_nuc();

    prob_t bg_pdf_edge_custom[NUC::SIZE];
    for(int j(0);j<NUC::SIZE;++j){
      bg_pdf_edge_custom[j] = r.bg_pdistro_nuc()[j]*(1.-edge_strength);
    }
    bg_pdf_edge_custom[n_edge] += edge_strength;

    if(is_edge_before_overedge){ bg_pdf_pos_before_overedge = bg_pdf_edge_custom; }
    else                       { bg_pdf_pos_after_overedge = bg_pdf_edge_custom; }

    prob_t bg_trans_rate[NUC::SIZE*NUC::SIZE];
    for(int n_overedge_1(0);n_overedge_1<NUC::SIZE;++n_overedge_1){
      for(int n_overedge_2(0);n_overedge_2<NUC::SIZE;++n_overedge_2){
        if(n_overedge_1==n_overedge_2) continue;
        smlfloat this_rate = r.get_category_nuc_context_mut_rate(NUC::N,static_cast<NUC::index_t>(n_overedge_1),
                                                                 NUC::N,static_cast<NUC::index_t>(n_overedge_2),
                                                                 site_category,group_category,
                                                                 bg_pdf_pos_before_overedge,
                                                                 bg_pdf_pos_after_overedge);

        bg_trans_rate[n_overedge_2+n_overedge_1*NUC::SIZE] = this_rate;
      }
    }
    matrix_balance_diagonal(bg_trans_rate,NUC::SIZE);
    prob_t bg_nuc_stat_pdf[NUC::SIZE];
    for(int j(0);j<NUC::SIZE;++j) { bg_nuc_stat_pdf[j] = r.bg_pdistro_nuc()[j]; }
    get_stationary_pdistro_from_rates(bg_nuc_stat_pdf,bg_trans_rate,NUC::SIZE,true);

    for(int j(0);j<NUC::SIZE;++j) {
      p_edge_overedge_tmp[n_edge][j] = bg_nuc_stat_pdf[j];
    }
  }

  finish_edge_dependency_distro(r.bg_pdistro_nuc(),
                                r.bg_pdistro_nuc(),
                                p_edge_overedge_tmp,
                                bg_pdf_nuc_overedge_cond_on_edge,
                                bg_pdf_nuc_edge_cond_on_overedge);
}




void
get_rate_edge_dependencies_noncoding_bidirectional_context(const rate_gtor_nuc_base& r,
                                                           const unsigned site_category,
                                                           const unsigned group_category,
                                                           const int overedge_dir,
                                                           prob_t bg_pdf_nuc_overedge_cond_on_edge_underedge[NUC::SIZE][NUC::SIZE][NUC::SIZE]){

  const bool is_edge_before_overedge(overedge_dir>0);

  smlfloat edge_tmp1,edge_tmp2,edge_tmp3;
  if(is_edge_before_overedge){
    edge_tmp1=r.edge_strength_pre();
    edge_tmp2=r.edge_strength_post();
    edge_tmp3=r.edge_strength_pre2();
  }//log_os << "PRE " << edge_tmp << "\n";}
  else                       {
    edge_tmp1=r.edge_strength_post();
    edge_tmp2=r.edge_strength_pre();
    edge_tmp3=r.edge_strength_post2();
  }//log_os << "POS " << edge_tmp << "\n";}

  const smlfloat edge_underedge_strength(edge_tmp3);
  const smlfloat edge_overedge_strength(edge_tmp1);
  const smlfloat edge_overedge_strength_back(edge_tmp2);

  // indexes correspond to: underedge,edge,overedge
  prob_t p_edge_tmp[NUC::SIZE][NUC::SIZE][NUC::SIZE];

  prob_t bg_pdf_edge_custom[NUC::SIZE][NUC::SIZE];
  for(unsigned i(0);i<NUC::SIZE;++i){
    for(unsigned j(0);j<NUC::SIZE;++j){
      bg_pdf_edge_custom[i][j] = r.bg_pdistro_nuc()[j]*(1.-edge_overedge_strength);
    }
    bg_pdf_edge_custom[i][i] += edge_overedge_strength;
  }

  prob_t bg_pdf_edge_custom_back[NUC::SIZE][NUC::SIZE];
  for(unsigned i(0);i<NUC::SIZE;++i){
    for(unsigned j(0);j<NUC::SIZE;++j){
      bg_pdf_edge_custom_back[i][j] = r.bg_pdistro_nuc()[j]*(1.-edge_overedge_strength_back);
    }
    bg_pdf_edge_custom_back[i][i] += edge_overedge_strength_back;
  }

  for(unsigned n_underedge(0);n_underedge<NUC::SIZE;++n_underedge){

    const prob_t* bg_pdf_pos_minus1 = r.bg_pdistro_nuc();
    const prob_t* bg_pdf_pos_plus1 = r.bg_pdistro_nuc();

    prob_t bg_pdf_underedge_custom[NUC::SIZE];
    for(int j(0);j<NUC::SIZE;++j){
      bg_pdf_underedge_custom[j] = r.bg_pdistro_nuc()[j]*(1.-edge_underedge_strength);
    }
    bg_pdf_underedge_custom[n_underedge] += edge_underedge_strength;

    if(is_edge_before_overedge){ bg_pdf_pos_minus1 = bg_pdf_underedge_custom; }
    else                       { bg_pdf_pos_plus1 = bg_pdf_underedge_custom; }

    prob_t bg_trans_rate[DINUC::SIZE*DINUC::SIZE];
    for(int d_eoe_1(0);d_eoe_1<DINUC::SIZE;++d_eoe_1){
      const NUC::index_t d1nx = DINUC::decode_nx(d_eoe_1);
      const NUC::index_t d1n0 = DINUC::decode_n0(d_eoe_1);

      for(int d_eoe_2(0);d_eoe_2<DINUC::SIZE;++d_eoe_2){
        if(d_eoe_1==d_eoe_2) continue;

        const NUC::index_t d2nx = DINUC::decode_nx(d_eoe_2);
        const NUC::index_t d2n0 = DINUC::decode_n0(d_eoe_2);

        unsigned nsubs(0);
        NUC::index_t d1ndiff(NUC::N);
        NUC::index_t d1ndiff_minus1(NUC::N);
        NUC::index_t d1ndiff_plus1(NUC::N);
        NUC::index_t d2ndiff(NUC::N);

        if( d1nx != d2nx ) {
          nsubs++; d1ndiff=d1nx; d2ndiff=d2nx;
          d1ndiff_plus1  = d1n0;
        }
        if( d1n0 != d2n0 ) {
          nsubs++; d1ndiff=d1n0; d2ndiff=d2n0;
          d1ndiff_minus1 = d1nx;
        }

        smlfloat this_rate(0.);
        if( nsubs == 1 ) {
          if(is_edge_before_overedge){
            if(d1nx != d2nx){
              bg_pdf_pos_minus1 = bg_pdf_underedge_custom;
              bg_pdf_pos_plus1 = bg_pdf_edge_custom[d1ndiff_plus1];
              d1ndiff_plus1 = NUC::N;
            } else {
#ifdef FOO13
              bg_pdf_pos_minus1 = r.bg_pdistro_nuc();
#else
              bg_pdf_pos_minus1 = bg_pdf_edge_custom_back[d1ndiff_minus1];
              d1ndiff_minus1 = NUC::N;
#endif
              bg_pdf_pos_plus1 = r.bg_pdistro_nuc();
            }
          } else {
            if( d1n0 != d2n0 ) {
              bg_pdf_pos_plus1 = bg_pdf_underedge_custom;
              bg_pdf_pos_minus1 = bg_pdf_edge_custom[d1ndiff_minus1];
              d1ndiff_minus1 = NUC::N;
            } else {
#ifdef FOO13
              bg_pdf_pos_plus1 = r.bg_pdistro_nuc();
#else
              bg_pdf_pos_plus1 = bg_pdf_edge_custom_back[d1ndiff_plus1];
              d1ndiff_plus1 = NUC::N;
#endif
              bg_pdf_pos_minus1 = r.bg_pdistro_nuc();
            }
          }

          this_rate = r.get_category_nuc_context_mut_rate(d1ndiff_minus1,d1ndiff,
                                                          d1ndiff_plus1,d2ndiff,
                                                          site_category,group_category,
                                                          bg_pdf_pos_minus1,
                                                          bg_pdf_pos_plus1);
        }
        bg_trans_rate[d_eoe_2+d_eoe_1*DINUC::SIZE] = this_rate;
      }
    }
    matrix_balance_diagonal(bg_trans_rate,DINUC::SIZE);
    prob_t bg_dinuc_stat_pdf[DINUC::SIZE];
    for(unsigned j(0);j<DINUC::SIZE;++j) {
      NUC::index_t nx = DINUC::decode_nx(j);
      NUC::index_t n0 = DINUC::decode_n0(j);
      bg_dinuc_stat_pdf[j] =  r.bg_pdistro_nuc()[nx]*r.bg_pdistro_nuc()[n0];
    }
    get_stationary_pdistro_from_rates(bg_dinuc_stat_pdf,bg_trans_rate,DINUC::SIZE,true);

    for(unsigned j(0);j<DINUC::SIZE;++j) {
      NUC::index_t nx = DINUC::decode_nx(j);
      NUC::index_t n0 = DINUC::decode_n0(j);
      if(is_edge_before_overedge){
        p_edge_tmp[n_underedge][nx][n0] = bg_dinuc_stat_pdf[j];
      } else {
        p_edge_tmp[n_underedge][n0][nx] = bg_dinuc_stat_pdf[j];
      }
    }
  }

  finish_edge_dependency_distro_bidirectional_context(r.bg_pdistro_nuc(),
                                                      r.bg_pdistro_nuc(),
                                                      p_edge_tmp,
                                                      bg_pdf_nuc_overedge_cond_on_edge_underedge);
}



#ifdef STAT_EDGE_CODING
void
get_rate_edge_dependencies_coding(const rate_gtor_nscodon_base& r,
                                  const rates_func_options& opt,
                                  const unsigned edge_pos,
                                  const int overedge_dir,
                                  prob_t* bg_pdf_nuc_overedge_cond_on_edge,
                                  prob_t* bg_pdf_nuc_edge_cond_on_overedge){

  const unsigned overedge_pos((edge_pos+(CODON::BASE_SIZE+overedge_dir))%CODON::BASE_SIZE);
  const bool is_edge_before_overedge(overedge_dir>0);

  const smlfloat edge_strength(is_edge_before_overedge ?
                               r.edge_strength_pre() :
                               r.edge_strength_post());

  const bool is_edge_in_overedge_codon(!((edge_pos==2 && is_edge_before_overedge) ||
                                         (edge_pos==0 && !is_edge_before_overedge)));

  prob_t p_edge_overedge_tmp[NUC::SIZE][NUC::SIZE];

#if 1 // use whole codon rate matrix to get stationary distro...

#if 0
  rates_func_options opt_nuc_pos_edge(opt);
  opt_nuc_pos_edge.is_use_custom_edge=true;
  opt_nuc_pos_edge.left_nuc_edge=r.bg_pdistro_nuc_pos(2);
  opt_nuc_pos_edge.right_nuc_edge=r.bg_pdistro_nuc_pos(0);
#endif
  die("nscodon stat based edge correction not (yet) working with experimental obs edge\n");

  smlfloat codon_rates[NSCODON::SIZE*NSCODON::SIZE];
  const prob_t* bg_pdistro_nuc_pos[CODON::BASE_SIZE];

  for(unsigned n_edge(0);n_edge<NUC::SIZE;++n_edge){

    prob_t bg_pdf_edge_custom[NUC::SIZE];
    for(int j(0);j<NUC::SIZE;++j){
      bg_pdf_edge_custom[j] = r.bg_pdistro_nuc_pos(edge_pos)[j]*(1.-edge_strength);
    }
    bg_pdf_edge_custom[n_edge] += edge_strength;

    for(unsigned i(0);i<CODON::BASE_SIZE;++i){
      bg_pdistro_nuc_pos[i] = r.bg_pdistro_nuc_pos(i);
    }

    if( ! is_edge_in_overedge_codon ){
      bg_pdistro_nuc_pos[edge_pos] = bg_pdf_edge_custom;
    }
#if 0
    rates_nscodon_context(r,SITE_MODEL::NSCODON,codon_rates,opt_nuc_pos_edge);
#endif
    prob_t codon_stat_pdf[NSCODON::SIZE];
    get_stationary_pdistro_from_rates(codon_stat_pdf,codon_rates,NSCODON::SIZE,false);


    for(unsigned j(0);j<NUC::SIZE;++j) p_edge_overedge_tmp[n_edge][j] = 0.;

    NUC::index_t jn[CODON::BASE_SIZE];
    for(unsigned j(0);j<NSCODON::SIZE;++j){
      SITE_MODEL::decode_nuc(SITE_MODEL::NSCODON,j,jn);
      smlfloat val( codon_stat_pdf[j]);
      if( is_edge_in_overedge_codon ){ val *= bg_pdf_edge_custom[jn[edge_pos]]; }
      p_edge_overedge_tmp[n_edge][jn[overedge_pos]] += val;
    }

    pdistro_norm(p_edge_overedge_tmp[n_edge],p_edge_overedge_tmp[n_edge]+NUC::SIZE);
  }


#else  // only use adjoining nucleotide stationary distro

  const unsigned overedge_codon_pos_before(periodic_decrement(overedge_pos,CODON::BASE_SIZE));
  const unsigned overedge_codon_pos_after(periodic_increment(overedge_pos,CODON::BASE_SIZE));

  for(unsigned n_edge(0);n_edge<NUC::SIZE;++n_edge){

    prob_t bg_pdf_edge_custom[NUC::SIZE];
    for(int j(0);j<NUC::SIZE;++j){
      bg_pdf_edge_custom[j] = r.bg_pdistro_nuc_pos(edge_pos)[j]*(1.-edge_strength);
    }
    bg_pdf_edge_custom[n_edge] += edge_strength;

    // setup plus1/minus1 distro's used to get the mutation rate:
    const prob_t* bg_pdf_pos_before_overedge(r.bg_pdistro_nuc_pos(overedge_codon_pos_before));
    const prob_t* bg_pdf_pos_after_overedge(r.bg_pdistro_nuc_pos(overedge_codon_pos_after));
    if(is_edge_before_overedge){ bg_pdf_pos_before_overedge = bg_pdf_edge_custom; }
    else                       { bg_pdf_pos_after_overedge = bg_pdf_edge_custom; }

    // setup codon position distro's used to calculate selection on a partial codon:
    const prob_t* bg_nuc_codon_pos[CODON::BASE_SIZE];
    for(unsigned i(0);i<CODON::BASE_SIZE;++i){ bg_nuc_codon_pos[i] = r.bg_pdistro_nuc_pos(i); }
    if(is_edge_in_overedge_codon){
      bg_nuc_codon_pos[edge_pos] = bg_pdf_edge_custom;
    }

    prob_t bg_trans_rate[NUC::SIZE*NUC::SIZE];
    for(int n_overedge_1(0);n_overedge_1<NUC::SIZE;++n_overedge_1){
      const NUC::index_t no1(static_cast<NUC::index_t>(n_overedge_1));

      for(int n_overedge_2(0);n_overedge_2<NUC::SIZE;++n_overedge_2){
        const NUC::index_t no2(static_cast<NUC::index_t>(n_overedge_2));
        if(no1 == no2) continue;

        smlfloat this_rate = r.get_category_nuc_context_mut_rate(NUC::N,no1,
                                                                 NUC::N,no2,
                                                                 opt.site_category,opt.group_category,
                                                                 bg_pdf_pos_before_overedge,
                                                                 bg_pdf_pos_after_overedge);

        // consider codon distro here for partial codon instead of nuc pos distros:
        this_rate *= get_average_codon_selection(r,no1,no2,overedge_pos,
                                                 bg_nuc_codon_pos,opt);

        bg_trans_rate[n_overedge_2+n_overedge_1*NUC::SIZE] = this_rate;
      }
    }
    matrix_balance_diagonal(bg_trans_rate,NUC::SIZE);
    prob_t bg_nuc_stat_pdf[NUC::SIZE];
    for(int j(0);j<NUC::SIZE;++j) { bg_nuc_stat_pdf[j] = r.bg_pdistro_nuc_pos(overedge_pos)[j]; }
    get_stationary_pdistro_from_rates(bg_nuc_stat_pdf,bg_trans_rate,NUC::SIZE,true);
    for(int j(0);j<NUC::SIZE;++j) {
      p_edge_overedge_tmp[n_edge][j] = bg_nuc_stat_pdf[j];
    }
  }
#endif


  finish_edge_dependency_distro(r.bg_pdistro_nuc_pos(edge_pos),
                                r.bg_pdistro_nuc_pos(overedge_pos),
                                p_edge_overedge_tmp,
                                bg_pdf_nuc_overedge_cond_on_edge,
                                bg_pdf_nuc_edge_cond_on_overedge);
}
#endif
#endif

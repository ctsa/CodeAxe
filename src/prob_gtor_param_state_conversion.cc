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

// $Id: prob_gtor_param_state_conversion.cc 1131 2008-01-29 04:08:32Z ctsa $

/// \file

#include "prob_gtor_param_state_conversion.h"
#include "stationary_pdistro.h"



prob_gtor_param_nsc4_seq_stat::
prob_gtor_param_nsc4_seq_stat(const bool is_c4_pre)
  : base_t(NSC4::SIZE,NSCODON::SIZE*NUC::SIZE), _is_c4_pre(is_c4_pre) {
  for(unsigned i(0);i<NUC::SIZE;++i){
    _is_free_block_start[NSCODON::SIZE*i] = true;
    _is_free_block_stop[NSCODON::SIZE*(i+1)-1] = true;
  }
}



void
prob_gtor_param_nsc4_seq_stat::
get_param_from_pdistro() {
  for(unsigned i(0);i<NUC::SIZE;++i){
    nsc4_pdf_2_conditional_nscodon_pdf(pdistro(),static_cast<NUC::index_t>(i),
                                       _param.begin()+NSCODON::SIZE*i);
  }
}



void
prob_gtor_param_nsc4_seq_stat::
param_norm() {
  for(unsigned i(0);i<NUC::SIZE;++i){
    pdistro_norm_free_param(_param.begin()+NSCODON::SIZE*i,
                            _param.begin()+NSCODON::SIZE*(i+1),
                            _is_train_param.begin()+NSCODON::SIZE*i);
  }
}



void
prob_gtor_param_nsc4_seq_stat::
get_pdistro_from_param() {

  prob_t nx_stat_distro[NUC::SIZE];

  // get stationary nx distro
  {
    const unsigned nx_codon_pos( _is_c4_pre ? 2 : 0 );

    // extract dependent nuc distro
    prob_t dep_nuc_p[NUC::SIZE*NUC::SIZE];
    for(unsigned i(0);i<NUC::SIZE;++i){
      nscodon_pdf_2_nuc_pdf_pos(_param.begin()+NSCODON::SIZE*i,
                                nx_codon_pos,dep_nuc_p+NUC::SIZE*i);
    }
    get_stationary_pdistro_from_pmatrix(nx_stat_distro,dep_nuc_p,NUC::SIZE);
  }

  // combine nx and conditional nuc distros:
  for(unsigned i(0);i<NSC4::SIZE;++i){
    const NSCODON::index_t c(NSC4::decode_nscodon(static_cast<NSC4::index_t>(i)));
    const NUC::index_t n(NSC4::decode_nuc(static_cast<NSC4::index_t>(i)));
    pdistro()[i] = nx_stat_distro[n]*_param[NSCODON::SIZE*n+c];
  }
}



/// adj_dinuc are the 2 bases closest to nx
///
template <typename RAI,typename RAI2>
void
nsc4_pdf_2_conditional_adj_dinuc_pdf(RAI nsc4_pdf,
                                     const NUC::index_t nx,
                                     const bool is_c4_pre,
                                     RAI2 dinuc_pdf){

  const unsigned adj_codon_pos( is_c4_pre ? 0 : 2 );
  NUC::index_t nc[CODON::BASE_SIZE];

  std::fill(dinuc_pdf,dinuc_pdf+DINUC::SIZE,0.);
  for(unsigned i(0);i<NSC4::SIZE;++i){
    const NSCODON::index_t c(NSC4::decode_nscodon(static_cast<NSC4::index_t>(i)));
    NSCODON::decode(nc,c);
    const NUC::index_t n(NSC4::decode_nuc(static_cast<NSC4::index_t>(i)));
    if(n!=nx) continue;

    const DINUC::index_t d(DINUC::encode(nc[1],nc[adj_codon_pos]));
    dinuc_pdf[d] += nsc4_pdf[i];
  }
  pdistro_norm(dinuc_pdf,dinuc_pdf+DINUC::SIZE);
}



prob_gtor_param_nsc4_seq_stat_nx2n::
prob_gtor_param_nsc4_seq_stat_nx2n(const bool is_c4_pre)
  : base_t(NSC4::SIZE,NSCODON::SIZE+TRINUC::SIZE), _is_c4_pre(is_c4_pre) {
  for(unsigned i(0);i<NUC::SIZE;++i){
    _is_free_block_start[DINUC::SIZE*i] = true;
    _is_free_block_stop[DINUC::SIZE*(i+1)-1] = true;
  }

  for(unsigned i(0);i<DINUC::SIZE;++i){
    const unsigned of(_dnc.get_offset(i,_is_c4_pre));
    _is_free_block_start[TRINUC::SIZE+of] = true;
    const unsigned si(_dnc.get_count(i,_is_c4_pre));
    _is_free_block_stop[TRINUC::SIZE+of+si-1] = true;
  }
}



void
prob_gtor_param_nsc4_seq_stat_nx2n::
get_param_from_pdistro() {
  for(unsigned i(0);i<NUC::SIZE;++i){
    nsc4_pdf_2_conditional_adj_dinuc_pdf(pdistro(),
                                         static_cast<NUC::index_t>(i),
                                         _is_c4_pre,
                                         _param.begin()+DINUC::SIZE*i);
  }

  prob_t nscodon_pdf[NSCODON::SIZE];
  nsc4_pdf_2_nscodon_pdf(pdistro(),nscodon_pdf);

  const unsigned adj_codon_pos( _is_c4_pre ? 0 : 2 );
  NUC::index_t nc[CODON::BASE_SIZE];

  prob_t pd[DINUC::SIZE];
  std::fill(pd,pd+DINUC::SIZE,0.);

  for(unsigned i(0);i<NSCODON::SIZE;++i){
    NSCODON::decode(nc,static_cast<NSCODON::index_t>(i));
    pd[ DINUC::encode(nc[1],nc[adj_codon_pos]) ] += nscodon_pdf[i];
  }

  for(unsigned i(0);i<NSCODON::SIZE;++i){
    const NSCODON::index_t c(static_cast<NSCODON::index_t>(i));
    NSCODON::decode(nc,c);
    const DINUC::index_t d(DINUC::encode(nc[1],nc[adj_codon_pos]));
    const unsigned in(_dnc.get_index(c,_is_c4_pre));
    _param[TRINUC::SIZE+in] = nscodon_pdf[i]/pd[d];
  }

#ifdef DEBUG
  static const smlfloat pthresh(1e-5);
  pdistro_check(pd,DINUC::SIZE,pthresh);

  for(unsigned i(0);i<DINUC::SIZE;++i){
    prob_t ptest(0.);

    // this loop checks sum of conditional codon distro, ie P(n1|n2,n3) for c4-post
    for(unsigned j(0);j<NSCODON::SIZE;++j){
      const NSCODON::index_t c(static_cast<NSCODON::index_t>(j));
      NSCODON::decode(nc,c);
      const DINUC::index_t d(DINUC::encode(nc[1],nc[adj_codon_pos]));
      if(d!=static_cast<DINUC::index_t>(i)) continue;
      const unsigned in(_dnc.get_index(c,_is_c4_pre));
      ptest += _param[TRINUC::SIZE+in];
    }
    assert(std::fabs(ptest-1.)<pthresh);

    // this loop checks the same thing in a different way:
    const unsigned of(_dnc.get_offset(i,_is_c4_pre));
    const unsigned si(_dnc.get_count(i,_is_c4_pre));
    pdistro_check(_param.begin()+TRINUC::SIZE+of,si,pthresh);
  }
#endif
}



void
prob_gtor_param_nsc4_seq_stat_nx2n::
param_norm() {
  for(unsigned i(0);i<NUC::SIZE;++i){
    pdistro_norm_free_param(_param.begin()+DINUC::SIZE*i,
                            _param.begin()+DINUC::SIZE*(i+1),
                            _is_train_param.begin()+DINUC::SIZE*i);
  }

  for(unsigned i(0);i<DINUC::SIZE;++i){
    const unsigned of(_dnc.get_offset(i,_is_c4_pre));
    const unsigned si(_dnc.get_count(i,_is_c4_pre));
    pdistro_norm_free_param(_param.begin()+TRINUC::SIZE+of,
                            _param.begin()+TRINUC::SIZE+of+si,
                            _is_train_param.begin()+TRINUC::SIZE+of);
  }
}



void
prob_gtor_param_nsc4_seq_stat_nx2n::
get_pdistro_from_param() {

  const unsigned adj_codon_pos( _is_c4_pre ? 0 : 2 );
  NUC::index_t nc[CODON::BASE_SIZE];

  prob_t nx_stat_distro[NUC::SIZE];

  // get stationary nx distro
  {
    const unsigned nx_codon_pos( _is_c4_pre ? 2 : 0 );

    // extract dependent nuc distro
    prob_t dep_nuc_p[NUC::SIZE*NUC::SIZE];
    for(unsigned i(0);i<NUC::SIZE;++i){
      const NUC::index_t nx = static_cast<NUC::index_t>(i);

      prob_t cond_nscodon_pdf[NSCODON::SIZE];
      for(unsigned j(0);j<NSCODON::SIZE;++j){
        const NSCODON::index_t c(static_cast<const NSCODON::index_t>(j));
        NSCODON::decode(nc,c);
        const DINUC::index_t d(DINUC::encode(nc[1],nc[adj_codon_pos]));
        const unsigned in(_dnc.get_index(c,_is_c4_pre));
        cond_nscodon_pdf[j] = _param[DINUC::SIZE*nx+d]*_param[TRINUC::SIZE+in];
      }

      nscodon_pdf_2_nuc_pdf_pos(cond_nscodon_pdf,
                                nx_codon_pos,dep_nuc_p+NUC::SIZE*i);
    }
    get_stationary_pdistro_from_pmatrix(nx_stat_distro,dep_nuc_p,NUC::SIZE);
  }

  // combine nx and conditional nuc distros:
  for(unsigned i(0);i<NSC4::SIZE;++i){
    const NSCODON::index_t c(NSC4::decode_nscodon(static_cast<NSC4::index_t>(i)));
    NSCODON::decode(nc,c);
    const NUC::index_t n(NSC4::decode_nuc(static_cast<NSC4::index_t>(i)));
    const DINUC::index_t d(DINUC::encode(nc[1],nc[adj_codon_pos]));
    const unsigned in(_dnc.get_index(c,_is_c4_pre));

    pdistro()[i] =
      nx_stat_distro[n]*_param[DINUC::SIZE*n+d]*_param[TRINUC::SIZE+in];
  }
}



// instantiate static class data:
//
const prob_gtor_param_nsc4_seq_stat_nx2n::dinuc_ns_count prob_gtor_param_nsc4_seq_stat_nx2n::_dnc;



prob_gtor_param_nsc4_seq_stat_nx2n::
dinuc_ns_count::
dinuc_ns_count(){
  NUC::index_t nc_pre[CODON::BASE_SIZE];
  NUC::index_t nc_post[CODON::BASE_SIZE];

  for(unsigned i(0);i<NUC::SIZE;++i){
    const NUC::index_t ni(static_cast<NUC::index_t>(i));
    nc_pre[0] = nc_post[2] = ni;

    for(unsigned j(0);j<NUC::SIZE;++j){
      const NUC::index_t nj(static_cast<NUC::index_t>(j));
      nc_pre[1] = nc_post[1] = nj;

      unsigned npre(0),npost(0);
      for(unsigned k(0);k<NUC::SIZE;++k){
        nc_pre[2] = nc_post[0] = static_cast<NUC::index_t>(k);

        if(NSCODON::encode(nc_pre) != NSCODON::NNN) npre++;
        if(NSCODON::encode(nc_post) != NSCODON::NNN) npost++;
      }

      const DINUC::index_t d(DINUC::encode(nj,ni));

      _count_pre[d] = npre;
      _count_post[d] = npost;
    }
  }

  _offset_pre[0] = 0;
  _offset_post[0] = 0;
  for(unsigned i(1);i<DINUC::SIZE;++i){
    _offset_pre[i] = _offset_pre[i-1]+_count_pre[i-1];
    _offset_post[i] = _offset_post[i-1]+_count_post[i-1];
  }

  // translate codon order:
  //
  for(unsigned i(0);i<NUC::SIZE;++i){
    const NUC::index_t ni(static_cast<NUC::index_t>(i));
    nc_pre[0] = nc_post[2] = ni;

    for(unsigned j(0);j<NUC::SIZE;++j){
      const NUC::index_t nj(static_cast<NUC::index_t>(j));
      nc_pre[1] = nc_post[1] = nj;

      const DINUC::index_t d(DINUC::encode(nj,ni));

      unsigned npre(0),npost(0);
      for(unsigned k(0);k<NUC::SIZE;++k){
        nc_pre[2] = nc_post[0] = static_cast<NUC::index_t>(k);

        const NSCODON::index_t cpre(NSCODON::encode(nc_pre));
        const NSCODON::index_t cpost(NSCODON::encode(nc_post));

        if(cpre != NSCODON::NNN) {
          _corder_pre[cpre] = get_offset(d,true)+npre;
          npre++;
        }
        if(cpost != NSCODON::NNN) {
         _corder_post[cpost] = get_offset(d,false)+npost;
          npost++;
        }
      }
    }
  }

#ifdef DEBUG
  unsigned tot1(0),tot2(0);
  for(unsigned i(0);i<DINUC::SIZE;++i){
    tot1 += get_count(i,true);
    tot2 += get_count(i,false);
  }
  assert(tot1 == NSCODON::SIZE && tot2 == NSCODON::SIZE);
#endif
}



/// adj nuc is the base closest to nx
///
template <typename RAI,typename RAI2>
void
nsc4_pdf_2_conditional_adj_nuc_pdf(RAI nsc4_pdf,
                                   const NUC::index_t nx,
                                   const bool is_c4_pre,
                                   RAI2 nuc_pdf){

  const unsigned adj_codon_pos( is_c4_pre ? 0 : 2 );
  NUC::index_t nc[CODON::BASE_SIZE];

  std::fill(nuc_pdf,nuc_pdf+NUC::SIZE,0.);
  for(unsigned i(0);i<NSC4::SIZE;++i){
    const NSCODON::index_t c(NSC4::decode_nscodon(static_cast<NSC4::index_t>(i)));
    NSCODON::decode(nc,c);
    const NUC::index_t n(NSC4::decode_nuc(static_cast<NSC4::index_t>(i)));
    if(n!=nx) continue;

    nuc_pdf[nc[adj_codon_pos]] += nsc4_pdf[i];
  }
  pdistro_norm(nuc_pdf,nuc_pdf+NUC::SIZE);
}



prob_gtor_param_nsc4_seq_stat_nx1n::
prob_gtor_param_nsc4_seq_stat_nx1n(const bool is_c4_pre)
  : base_t(NSC4::SIZE,NSCODON::SIZE+DINUC::SIZE), _is_c4_pre(is_c4_pre) {
  for(unsigned i(0);i<NUC::SIZE;++i){
    _is_free_block_start[NUC::SIZE*i] = true;
    _is_free_block_stop[NUC::SIZE*(i+1)-1] = true;
  }

  for(unsigned i(0);i<NUC::SIZE;++i){
    const NUC::index_t ni(static_cast<NUC::index_t>(i));
    const unsigned of(_nnc.get_offset(ni,_is_c4_pre));
    _is_free_block_start[DINUC::SIZE+of] = true;
    const unsigned si(_nnc.get_count(ni,_is_c4_pre));
    _is_free_block_stop[DINUC::SIZE+of+si-1] = true;
  }
}



void
prob_gtor_param_nsc4_seq_stat_nx1n::
get_param_from_pdistro() {
  for(unsigned i(0);i<NUC::SIZE;++i){
    nsc4_pdf_2_conditional_adj_nuc_pdf(pdistro(),
                                       static_cast<NUC::index_t>(i),
                                       _is_c4_pre,
                                       _param.begin()+NUC::SIZE*i);
  }

  prob_t nscodon_pdf[NSCODON::SIZE];
  nsc4_pdf_2_nscodon_pdf(pdistro(),nscodon_pdf);

  const unsigned adj_codon_pos( _is_c4_pre ? 0 : 2 );
  NUC::index_t nc[CODON::BASE_SIZE];

  prob_t pn[NUC::SIZE];
  std::fill(pn,pn+NUC::SIZE,0.);

  for(unsigned i(0);i<NSCODON::SIZE;++i){
    NSCODON::decode(nc,static_cast<NSCODON::index_t>(i));
    pn[ nc[adj_codon_pos] ] += nscodon_pdf[i];
  }

  for(unsigned i(0);i<NSCODON::SIZE;++i){
    const NSCODON::index_t c(static_cast<NSCODON::index_t>(i));
    NSCODON::decode(nc,c);
    const unsigned in(_nnc.get_index(c,_is_c4_pre));
    _param[DINUC::SIZE+in] = nscodon_pdf[i]/pn[nc[adj_codon_pos]];
  }

#ifdef DEBUG
  static const smlfloat pthresh(1e-5);
  pdistro_check(pn,NUC::SIZE,pthresh);

  for(unsigned i(0);i<NUC::SIZE;++i){
    const NUC::index_t ni(static_cast<NUC::index_t>(i));

    prob_t ptest(0.);
    // this loop checks the sum of the conditional codon distro, ie P(n1,n2|n3) for c4-post
    for(unsigned j(0);j<NSCODON::SIZE;++j){
      const NSCODON::index_t c(static_cast<NSCODON::index_t>(j));
      NSCODON::decode(nc,c);
      if(nc[adj_codon_pos]!=ni) continue;
      const unsigned in(_nnc.get_index(c,_is_c4_pre));
      ptest += _param[DINUC::SIZE+in];
    }
    assert(std::fabs(ptest-1.)<pthresh);

    // this checks the same thing using a different method
    const unsigned of(_nnc.get_offset(ni,_is_c4_pre));
    const unsigned si(_nnc.get_count(ni,_is_c4_pre));
    pdistro_check(_param.begin()+DINUC::SIZE+of,si,pthresh);
  }
#endif
}



void
prob_gtor_param_nsc4_seq_stat_nx1n::
param_norm() {
  for(unsigned i(0);i<NUC::SIZE;++i){
    pdistro_norm_free_param(_param.begin()+NUC::SIZE*i,
                            _param.begin()+NUC::SIZE*(i+1),
                            _is_train_param.begin()+NUC::SIZE*i);
  }

  for(unsigned i(0);i<NUC::SIZE;++i){
    const NUC::index_t ni(static_cast<NUC::index_t>(i));

    const unsigned of(_nnc.get_offset(ni,_is_c4_pre));
    const unsigned si(_nnc.get_count(ni,_is_c4_pre));
    pdistro_norm_free_param(_param.begin()+DINUC::SIZE+of,
                            _param.begin()+DINUC::SIZE+of+si,
                            _is_train_param.begin()+DINUC::SIZE+of);
  }
}



void
prob_gtor_param_nsc4_seq_stat_nx1n::
get_pdistro_from_param() {

  const unsigned adj_codon_pos( _is_c4_pre ? 0 : 2 );
  NUC::index_t nc[CODON::BASE_SIZE];

  prob_t nx_stat_distro[NUC::SIZE];

  // get stationary nx distro
  {
    const unsigned nx_codon_pos( _is_c4_pre ? 2 : 0 );

    // extract dependent nuc distro
    prob_t dep_nuc_p[NUC::SIZE*NUC::SIZE];
    for(unsigned i(0);i<NUC::SIZE;++i){
      const NUC::index_t nx = static_cast<NUC::index_t>(i);

      prob_t cond_nscodon_pdf[NSCODON::SIZE];
      for(unsigned j(0);j<NSCODON::SIZE;++j){
        const NSCODON::index_t c(static_cast<const NSCODON::index_t>(j));
        NSCODON::decode(nc,c);
        const unsigned in(_nnc.get_index(c,_is_c4_pre));
        cond_nscodon_pdf[j] = _param[NUC::SIZE*nx+nc[adj_codon_pos]]*_param[DINUC::SIZE+in];
      }

      nscodon_pdf_2_nuc_pdf_pos(cond_nscodon_pdf,
                                nx_codon_pos,dep_nuc_p+NUC::SIZE*i);
    }
    get_stationary_pdistro_from_pmatrix(nx_stat_distro,dep_nuc_p,NUC::SIZE);
  }

  // combine nx and conditional nuc distros:
  for(unsigned i(0);i<NSC4::SIZE;++i){
    const NSCODON::index_t c(NSC4::decode_nscodon(static_cast<NSC4::index_t>(i)));
    NSCODON::decode(nc,c);
    const NUC::index_t n(NSC4::decode_nuc(static_cast<NSC4::index_t>(i)));
    const unsigned in(_nnc.get_index(c,_is_c4_pre));

    pdistro()[i] =
      nx_stat_distro[n]*_param[NUC::SIZE*n+nc[adj_codon_pos]]*_param[DINUC::SIZE+in];
  }
}



// instantiate static class data:
//
const prob_gtor_param_nsc4_seq_stat_nx1n::nuc_ns_count prob_gtor_param_nsc4_seq_stat_nx1n::_nnc;



prob_gtor_param_nsc4_seq_stat_nx1n::
nuc_ns_count::
nuc_ns_count(){
  NUC::index_t nc_pre[CODON::BASE_SIZE];
  NUC::index_t nc_post[CODON::BASE_SIZE];

  for(unsigned i(0);i<NUC::SIZE;++i){
    const NUC::index_t ni(static_cast<NUC::index_t>(i));
    nc_pre[0] = nc_post[2] = ni;

    unsigned npre(0),npost(0);
    for(unsigned j(0);j<NUC::SIZE;++j){
      const NUC::index_t nj(static_cast<NUC::index_t>(j));
      nc_pre[1] = nc_post[1] = nj;

      for(unsigned k(0);k<NUC::SIZE;++k){
        nc_pre[2] = nc_post[0] = static_cast<NUC::index_t>(k);

        if(NSCODON::encode(nc_pre) != NSCODON::NNN) npre++;
        if(NSCODON::encode(nc_post) != NSCODON::NNN) npost++;
      }
    }
    _count_pre[ni] = npre;
    _count_post[ni] = npost;
  }

  _offset_pre[0] = 0;
  _offset_post[0] = 0;
  for(unsigned i(1);i<NUC::SIZE;++i){
    _offset_pre[i] = _offset_pre[i-1]+_count_pre[i-1];
    _offset_post[i] = _offset_post[i-1]+_count_post[i-1];
  }

  // translate codon order:
  //
  for(unsigned i(0);i<NUC::SIZE;++i){
    const NUC::index_t ni(static_cast<NUC::index_t>(i));
    nc_pre[0] = nc_post[2] = ni;

    unsigned npre(0),npost(0);
    for(unsigned j(0);j<NUC::SIZE;++j){
      nc_pre[1] = nc_post[1] = static_cast<NUC::index_t>(j);

      for(unsigned k(0);k<NUC::SIZE;++k){
        nc_pre[2] = nc_post[0] = static_cast<NUC::index_t>(k);

        const NSCODON::index_t cpre(NSCODON::encode(nc_pre));
        const NSCODON::index_t cpost(NSCODON::encode(nc_post));

        if(cpre != NSCODON::NNN) {
          _corder_pre[cpre] = get_offset(ni,true)+npre;
          npre++;
        }
        if(cpost != NSCODON::NNN) {
         _corder_post[cpost] = get_offset(ni,false)+npost;
          npost++;
        }
      }
    }
  }

#ifdef DEBUG
  unsigned tot1(0),tot2(0);
  for(unsigned i(0);i<NUC::SIZE;++i){
    const NUC::index_t ni(static_cast<NUC::index_t>(i));

    tot1 += get_count(ni,true);
    tot2 += get_count(ni,false);
  }
  assert(tot1 == NSCODON::SIZE && tot2 == NSCODON::SIZE);
#endif
}

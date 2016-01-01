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

// $Id: prob_gtor_param_state_conversion.h 1083 2008-01-07 22:41:26Z ctsa $

/// \file

#ifndef __PROB_GTOR_PARAM_STATE_CONVERSION_H
#define __PROB_GTOR_PARAM_STATE_CONVERSION_H

#include "prob_gtor_param.h"
#include "util/bio/bioseq_util.h"
#include "util/bio/bioseq_util_pdf_conversions.h"



#if 0
/// generalization of the distros below -- in progress
///
struct prob_gtor_param_bio_pdistro_builder : public prob_gtor_param {

  typedef prob_param_gtor base_t;
  typedef prob_gtor_bio_pdistro_builder self_t;

  prob_gtor_bio_pdistro_builder(const SITE_MODEL::index_t target,
                                const BIO_PDISTRO::index_t param)
    : base_t(BIO_PDISTRO::state_size(BIO_PDISTRO::convert_from_site_model(target)),
             BIO_PDISTRO::state_size(param)) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void get_param_from_pdistro() {
    bio_pdistro_convert(BIO_PDISTRO::convert_from_site_model(target),param,pdistro(),x);
  }

  virtual
  void pdistro_norm() {
    const unsigned ps(param_size());
    FloatType* x2(new FloatType[ps]);
    bio_pdistro_norm(param,x,x2);
  }

  virtual
  void get_pdistro_from_param() {
    const unsigned ps(param_size());
    FloatType* x2(new FloatType[ps]);
    bio_pdistro_convert(param,BIO_PDISTRO::convert_from_site_model(target),x2,pdistro());
    delete [] x2;
  }
};
#endif



struct prob_gtor_param_dinuc_from_nuc : public prob_gtor_param {

  typedef prob_gtor_param base_t;
  typedef prob_gtor_param_dinuc_from_nuc self_t;

  prob_gtor_param_dinuc_from_nuc() : base_t(DINUC::SIZE,NUC::SIZE) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void get_param_from_pdistro() {
    dinuc_pdf_2_nuc_pdf_pos1(pdistro(),_param.begin());
  }

  virtual
  void get_pdistro_from_param() {
    nuc_pdf_2_dinuc_pdf(_param.begin(),pdistro());
  }
};



/// \brief an nsc4 distro enforcing sequence stationarity, so that the
/// marginal distro of (for nsc4-post, for example) n1 and n4 must be
/// equal
///
struct prob_gtor_param_nsc4_seq_stat : public prob_gtor_param {

  typedef prob_gtor_param base_t;
  typedef prob_gtor_param_nsc4_seq_stat self_t;

  explicit
  prob_gtor_param_nsc4_seq_stat(const bool is_c4_pre);

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual void get_param_from_pdistro();

  virtual void param_norm();

  virtual void get_pdistro_from_param();

  bool _is_c4_pre;
};



/// an nsc4 distro enforcing sequence stationarity, so that the
/// marginal distro of (for nsc4-post, for example) n1 and n4 must be
/// equal, this one does not allow any direct dependency between n1
/// and n4, instead using P(n2,n3|n4)*P(n1|n2,n3) to get the
/// conditional codon distro, from which the c4 disto is found by the
/// stationary marginal distro for n1/n4
///
struct prob_gtor_param_nsc4_seq_stat_nx2n : public prob_gtor_param {

  typedef prob_gtor_param base_t;
  typedef prob_gtor_param_nsc4_seq_stat_nx2n self_t;

  explicit
  prob_gtor_param_nsc4_seq_stat_nx2n(const bool is_c4_pre);

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual void get_param_from_pdistro();

  virtual void param_norm();

  virtual void get_pdistro_from_param();

  struct dinuc_ns_count {

    dinuc_ns_count();

    unsigned get_count(const DINUC::index_t d,
                       const bool is_c4_pre) const {
      if(is_c4_pre) return _count_pre[d];
      else          return _count_post[d];
    }

    unsigned get_offset(const DINUC::index_t d,
                        const bool is_c4_pre) const {
      if(is_c4_pre) return _offset_pre[d];
      else          return _offset_post[d];
    }

    unsigned get_index(const NSCODON::index_t c,
                       const bool is_c4_pre) const {
      if(is_c4_pre) return _corder_pre[c];
      else          return _corder_post[c];
    }

  private:
    unsigned _count_pre[DINUC::SIZE];
    unsigned _count_post[DINUC::SIZE];
    unsigned _offset_pre[DINUC::SIZE];
    unsigned _offset_post[DINUC::SIZE];
    unsigned _corder_pre[NSCODON::SIZE];
    unsigned _corder_post[NSCODON::SIZE];
  };

  static const dinuc_ns_count _dnc;
  bool _is_c4_pre;
};


/// an nsc4 distro enforcing sequence stationarity, so that the
/// marginal distro of (for nsc4-post, for example) n1 and n4 must be
/// equal, this one does not allow any direct dependency between
/// (n1,n2) and n4, instead using P(n3|n4)*P(n1,n2|n3) to get the
/// conditional codon distro, from which the c4 disto is found by the
/// stationary marginal distro for n1/n4
///
struct prob_gtor_param_nsc4_seq_stat_nx1n : public prob_gtor_param {

  typedef prob_gtor_param base_t;
  typedef prob_gtor_param_nsc4_seq_stat_nx1n self_t;

  explicit
  prob_gtor_param_nsc4_seq_stat_nx1n(const bool is_c4_pre);

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual void get_param_from_pdistro();

  virtual void param_norm();

  virtual void get_pdistro_from_param();

  struct nuc_ns_count {

    nuc_ns_count();

    unsigned get_count(const NUC::index_t n,
                       const bool is_c4_pre) const {
      if(is_c4_pre) return _count_pre[n];
      else          return _count_post[n];
    }

    unsigned get_offset(const NUC::index_t n,
                        const bool is_c4_pre) const {
      if(is_c4_pre) return _offset_pre[n];
      else          return _offset_post[n];
    }

    unsigned get_index(const NSCODON::index_t c,
                       const bool is_c4_pre) const {
      if(is_c4_pre) return _corder_pre[c];
      else          return _corder_post[c];
    }

  private:
    unsigned _count_pre[NUC::SIZE];
    unsigned _count_post[NUC::SIZE];
    unsigned _offset_pre[NUC::SIZE];
    unsigned _offset_post[NUC::SIZE];
    unsigned _corder_pre[NSCODON::SIZE];
    unsigned _corder_post[NSCODON::SIZE];
  };

  static const nuc_ns_count _nnc;
  bool _is_c4_pre;
};



/// \brief nsc4 distro built from nscodon parameters
///
struct prob_gtor_param_nsc4_from_nscodon : public prob_gtor_param {

  typedef prob_gtor_param base_t;
  typedef prob_gtor_param_nsc4_from_nscodon self_t;

  explicit
  prob_gtor_param_nsc4_from_nscodon(const bool is_c4_pre)
    : base_t(NSC4::SIZE,NSCODON::SIZE), _is_c4_pre(is_c4_pre) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void get_param_from_pdistro() {
    nsc4_pdf_2_nscodon_pdf(pdistro(),_param.begin());
  }

  virtual
  void get_pdistro_from_param() {
    if(_is_c4_pre){
      nscodon_pdf_2_nsc4_pre_pdf(_param.begin(),pdistro());
    } else {
      nscodon_pdf_2_nsc4_post_pdf(_param.begin(),pdistro());
    }
  }

  bool _is_c4_pre;
};




/// \brief nsc4 distro build from positional nuc params
///
struct prob_gtor_param_nsc4_from_nuc_pos : public prob_gtor_param {

  typedef prob_gtor_param base_t;
  typedef prob_gtor_param_nsc4_from_nuc_pos self_t;

  explicit
  prob_gtor_param_nsc4_from_nuc_pos(const bool is_c4_pre)
    : base_t(NSC4::SIZE,PARAM_SIZE), _is_c4_pre(is_c4_pre) {
    for(unsigned i(1);i<CODON::BASE_SIZE;++i){
      _is_free_block_stop[i*NUC::SIZE-1] = true;
      _is_free_block_start[i*NUC::SIZE] = true;
    }
  }

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void get_param_from_pdistro() {
    nsc4_pdf_2_nuc_pdf_all_pos(pdistro(),_param.begin());
  }

  virtual
  void param_norm() {
    for(unsigned i(0);i<CODON::BASE_SIZE;++i){
      pdistro_norm_free_param(_param.begin()+NUC::SIZE*i,
                              _param.begin()+NUC::SIZE*(i+1),
                              _is_train_param.begin()+NUC::SIZE*i);
    }
  }

  virtual
  void get_pdistro_from_param() {
    if(_is_c4_pre){
      nuc_pdf_all_pos_2_nsc4_pre_pdf(_param.begin(),pdistro());
    } else {
      nuc_pdf_all_pos_2_nsc4_post_pdf(_param.begin(),pdistro());
    }
  }

  enum { PARAM_SIZE = NUC::SIZE*CODON::BASE_SIZE };

  bool _is_c4_pre;
};




/// \brief nsc4 distro (pre or post) build from nuc params
///
struct prob_gtor_param_nsc4_from_nuc : public prob_gtor_param {

  typedef prob_gtor_param base_t;
  typedef prob_gtor_param_nsc4_from_nuc self_t;

  prob_gtor_param_nsc4_from_nuc() : base_t(NSC4::SIZE,NUC::SIZE) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void get_param_from_pdistro() {
    nsc4_pdf_2_nuc_pdf(pdistro(),_param.begin());
  }

  virtual
  void get_pdistro_from_param() {
    nuc_pdf_2_nsc4_pdf(_param.begin(),pdistro());
  }
};



/// \brief frequencies for nsc5 states build from
/// nscodon frequency parameters
///
struct prob_gtor_param_nsc5_from_nscodon : public prob_gtor_param {

  typedef prob_gtor_param base_t;
  typedef prob_gtor_param_nsc5_from_nscodon self_t;

  prob_gtor_param_nsc5_from_nscodon() : base_t(NSC5::SIZE,NSCODON::SIZE) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void get_param_from_pdistro() {
    nsc5_pdf_2_nscodon_pdf(pdistro(),_param.begin());
  }

  virtual
  void get_pdistro_from_param() {
    nscodon_pdf_2_nsc5_pdf(_param.begin(),pdistro());
  }
};




/// \brief frequencies for nsc4 states (pre or post) build from
/// nscodon frequency parameters
///
struct prob_gtor_param_nscodon_from_nuc : public prob_gtor_param {

  typedef prob_gtor_param base_t;
  typedef prob_gtor_param_nscodon_from_nuc self_t;

  prob_gtor_param_nscodon_from_nuc() : base_t(NSCODON::SIZE,NUC::SIZE) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void get_param_from_pdistro() {
    nscodon_pdf_2_nuc_pdf(pdistro(),_param.begin());
  }

  virtual
  void get_pdistro_from_param() {
    nuc_pdf_2_nscodon_pdf(_param.begin(),pdistro());
  }
};




/// \brief frequencies for nsc4 states (pre or post) build from
/// nscodon frequency parameters
///
struct prob_gtor_param_nscodon_from_nuc_pos : public prob_gtor_param {

  typedef prob_gtor_param base_t;
  typedef prob_gtor_param_nscodon_from_nuc_pos self_t;

  prob_gtor_param_nscodon_from_nuc_pos() : base_t(NSCODON::SIZE,PARAM_SIZE) {
    for(unsigned i(1);i<CODON::BASE_SIZE;++i){
      _is_free_block_stop[i*NUC::SIZE-1] = true;
      _is_free_block_start[i*NUC::SIZE] = true;
    }
  }

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void get_param_from_pdistro() {
    nscodon_pdf_2_nuc_pdf_all_pos(pdistro(),_param.begin());
  }

  virtual
  void param_norm() {
    for(unsigned i(0);i<CODON::BASE_SIZE;++i){
      pdistro_norm_free_param(_param.begin()+NUC::SIZE*i,
                              _param.begin()+NUC::SIZE*(i+1),
                              _is_train_param.begin()+NUC::SIZE*i);
    }
  }

  virtual
  void get_pdistro_from_param() {
    nuc_pdf_all_pos_2_nscodon_pdf(_param.begin(),pdistro());
  }

  enum { PARAM_SIZE = NUC::SIZE*CODON::BASE_SIZE };
};



struct prob_gtor_param_trinuc_from_nuc : public prob_gtor_param {

  typedef prob_gtor_param base_t;

  typedef prob_gtor_param_trinuc_from_nuc self_t;

  prob_gtor_param_trinuc_from_nuc() : base_t(TRINUC::SIZE,NUC::SIZE) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void get_param_from_pdistro() {
    trinuc_pdf_2_nuc_pdf_pos(pdistro(),2,_param.begin());
  }

  virtual
  void get_pdistro_from_param() {
    nuc_pdf_2_trinuc_pdf(_param.begin(),pdistro());
  }
};


#endif

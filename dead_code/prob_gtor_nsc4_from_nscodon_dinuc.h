// -*- mode: c++; indent-tabs-mode: nil; -*-

/// \file prob_gtor_nsc4_from_nscodon_dinuc.h
/// \brief
///

#ifndef __PROB_GTOR_NSC4_FROM_NSCODON_DINUC_H
#define __PROB_GTOR_NSC4_FROM_NSCODON_DINUC_H


#include "bioseq_util.h"
#include "prob_gtor.h"
#include "rate_edge_dependencies.h"


/// \brief generate nsc4 distro from nscodon distro plus the n2,n0
/// codon boundary joint nuc distro
///
template <typename FLOAT_T>
struct prob_gtor_nsc4_from_nscodon_dinuc : public prob_gtor<FLOAT_T> {

  typedef prob_gtor<FLOAT_T> base_t;
  typedef prob_gtor_nsc4_from_nscodon_dinuc<FLOAT_T> self_t;

  prob_gtor_nsc4_from_nscodon_dinuc(const bool is_c4_pre)
    : base_t(NSC4::SIZE,NSCODON::SIZE+DINUC::SIZE), _is_c4_pre(is_c4_pre) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void param_state_internal(FLOAT_T* x) const {
    nsc4_pdf_2_nscodon_pdf(this->distro(),x);
    nsc4_pdf_2_codon_boundary_dinuc_pdf(this->distro(),_is_c4_pre,x+NSCODON::SIZE);

    shape_dinuc_nscodon(x,x+NSCODON::SIZE);
  }

  virtual
  void set_param_state_internal(const FLOAT_T* x) {
    FLOAT_T x_nscodon[NSCODON::SIZE];
    for(unsigned i(0);i<NSCODON::SIZE;++i) x_nscodon[i] = x[i];
    pdistro_norm(x_nscodon,NSCODON::SIZE);

    FLOAT_T x_dinuc[DINUC::SIZE];
    for(unsigned i(0);i<DINUC::SIZE;++i) x_dinuc[i] = x[i+NSCODON::SIZE];
    pdistro_norm(x_dinuc,DINUC::SIZE);

    shape_dinuc_nscodon(x_nscodon,x_dinuc);

    if(_is_c4_pre){
      nscodon_pdf_dinuc_pdf_2_nsc4_pre_pdf(x_nscodon,x_dinuc,this->distro());
    } else {
      nscodon_pdf_dinuc_pdf_2_nsc4_post_pdf(x_nscodon,x_dinuc,this->distro());
    }
  }

  // dinuc distro is used on which side of the codon:
  bool _is_c4_pre;
};


#endif

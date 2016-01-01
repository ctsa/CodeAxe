// -*- mode: c++; indent-tabs-mode: nil; -*-

/// \file prob_gtor_dinuc_from_nuc.h
/// \brief
///


#ifndef __PROB_GTOR_DINUC_FROM_NUC_H
#define __PROB_GTOR_DINUC_FROM_NUC_H

#include "bioseq_util.h"
#include "bioseq_util_pdf_conversions.h"
#include "prob_gtor.h"


/// \brief
///
template <typename FLOAT_T>
struct prob_gtor_dinuc_from_nuc : public prob_gtor<FLOAT_T> {

  typedef prob_gtor<FLOAT_T> base_t;
  typedef prob_gtor_dinuc_from_nuc<FLOAT_T> self_t;

  prob_gtor_dinuc_from_nuc() : base_t(DINUC::SIZE,PARAM_SIZE) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void param_state_internal(FLOAT_T* x) const {
    dinuc_pdf_2_nuc_pdf_pos1(this->distro(),x);
  }

  virtual
  void set_param_state_internal(const FLOAT_T* x) {
    FLOAT_T x2[PARAM_SIZE];
    for(unsigned i(0);i<PARAM_SIZE;++i) x2[i] = x[i];
    pdistro_norm(x2,PARAM_SIZE);
    nuc_pdf_2_dinuc_pdf(x2,this->distro());
  }

  enum { PARAM_SIZE = NUC::SIZE };
};


#endif

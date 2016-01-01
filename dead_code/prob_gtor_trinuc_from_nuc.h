// -*- mode: c++; indent-tabs-mode: nil; -*-

/// \file prob_gtor_trinuc_from_nuc.h
/// \brief
///

#ifndef __PROB_GTOR_TRINUC_FROM_NUC_H
#define __PROB_GTOR_TRINUC_FROM_NUC_H


#include "bioseq_util.h"
#include "bioseq_util_pdf_conversions.h"
#include "prob_gtor.h"


/// \brief frequencies for nsc4 states (pre or post) build from
/// nscodon frequency parameters
///
template <typename FLOAT_T>
struct prob_gtor_trinuc_from_nuc : public prob_gtor<FLOAT_T> {

  typedef prob_gtor<FLOAT_T> base_t;

  typedef prob_gtor_trinuc_from_nuc<FLOAT_T> self_t;

  prob_gtor_trinuc_from_nuc() : base_t(TRINUC::SIZE,PARAM_SIZE) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void param_state_internal(FLOAT_T* x) const {
    trinuc_pdf_2_nuc_pdf_pos(this->distro(),2,x);
  }

  virtual
  void set_param_state_internal(const FLOAT_T* x) {
    FLOAT_T x2[PARAM_SIZE];
    for(unsigned i(0);i<PARAM_SIZE;++i) x2[i] = x[i];
    pdistro_norm(x2,PARAM_SIZE);
    nuc_pdf_2_trinuc_pdf(x2,this->distro());
  }

  enum { PARAM_SIZE = NUC::SIZE };
};


#endif

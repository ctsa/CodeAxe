// -*- mode: c++; indent-tabs-mode: nil; -*-

/// \file prob_gtor_nsc5_from_nscodon.h
/// \brief
///

#ifndef __PROB_GTOR_NSC5_FROM_NSCODON_H
#define __PROB_GTOR_NSC5_FROM_NSCODON_H


#include "bioseq_util.h"
#include "bioseq_util_pdf_conversions.h"
#include "prob_gtor.h"


/// \brief frequencies for nsc5 states build from
/// nscodon frequency parameters
///
template <typename FLOAT_T>
struct prob_gtor_nsc5_from_nscodon : public prob_gtor<FLOAT_T> {

  typedef prob_gtor<FLOAT_T> base_t;
  typedef prob_gtor_nsc5_from_nscodon<FLOAT_T> self_t;

  prob_gtor_nsc5_from_nscodon() : base_t(NSC5::SIZE,PARAM_SIZE) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void param_state_internal(FLOAT_T* x) const {
    nsc5_pdf_2_nscodon_pdf(this->distro(),x);
  }

  virtual
  void set_param_state_internal(const FLOAT_T* x) {
    FLOAT_T x2[PARAM_SIZE];
    for(unsigned i(0);i<PARAM_SIZE;++i) x2[i] = x[i];
    pdistro_norm(x2,PARAM_SIZE);
    nscodon_pdf_2_nsc5_pdf(x2,this->distro());
  }

  enum { PARAM_SIZE = NSCODON::SIZE };
};



#endif

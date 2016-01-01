// -*- mode: c++; indent-tabs-mode: nil; -*-

/// \file prob_gtor_nscodon_from_nuc_pos.h
/// \brief
///

#ifndef __PROB_GTOR_NSCODON_FROM_NUC_POS_H
#define __PROB_GTOR_NSCODON_FROM_NUC_POS_H


#include "bioseq_util.h"
#include "bioseq_util_pdf_conversions.h"
#include "prob_gtor.h"


/// \brief frequencies for nsc4 states (pre or post) build from
/// nscodon frequency parameters
///
template <typename FLOAT_T>
struct prob_gtor_nscodon_from_nuc_pos : public prob_gtor<FLOAT_T> {

  typedef prob_gtor<FLOAT_T> base_t;
  typedef prob_gtor_nscodon_from_nuc_pos<FLOAT_T> self_t;

  prob_gtor_nscodon_from_nuc_pos() : base_t(NSCODON::SIZE,PARAM_SIZE) {}

  virtual base_t* clone() const { return new self_t(*this); }

private:
  virtual
  void param_state_internal(FLOAT_T* x) const {
    nscodon_pdf_2_nuc_pdf_all_pos(this->distro(),x);
  }

  virtual
  void set_param_state_internal(const FLOAT_T* x) {
    FLOAT_T x2[PARAM_SIZE];
    for(unsigned i(0);i<PARAM_SIZE;++i) x2[i] = x[i];
    for(unsigned i(0);i<CODON::BASE_SIZE;++i){
      pdistro_norm(x2+NUC::SIZE*i,NUC::SIZE);
    }
    nuc_pdf_all_pos_2_nscodon_pdf(x2,this->distro());
  }

  enum { PARAM_SIZE = NUC::SIZE*CODON::BASE_SIZE };
};


#endif

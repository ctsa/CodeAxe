// -*- mode: c++; indent-tabs-mode: nil; -*-

/// \file rate_gtor_nsc4.cc
/// \brief
///


#include "rate_gtor_nsc4.h"
#include "util/bioseq_util_flux_rate_conversions.h"
#include "util/die.h"


rate_gtor_nsc4::
rate_gtor_nsc4(const rate_gtor_nuc_options nopt,
               const rate_gtor_nscodon_options copt,
               const RATE_GTOR_MODEL::index_t c4_type)
  : base_t(nopt,copt,NSC4::SIZE), _c4_type(c4_type) {
  if( _c4_type != RATE_GTOR_MODEL::C4PRE && _c4_type != RATE_GTOR_MODEL::C4POST ) {
    die("Invalid c4 type for nsc4 rate_gtor model");
  }
  bg_pdistro_update();
}


void
rate_gtor_nsc4::
flux_nscodon_no_aa_selection(smlfloat nscodon_flux[NSCODON::SIZE*NSCODON::SIZE]) const {

  smlfloat* nsc4_rates(new smlfloat[NSC4::SIZE*NSC4::SIZE]);
  rates(nsc4_rates,rates_func_options(-1,-1,true));
  nsc4_rates_2_nscodon_flux(nscodon_flux,nsc4_rates,bg_pdistro());
  delete [] nsc4_rates;
}

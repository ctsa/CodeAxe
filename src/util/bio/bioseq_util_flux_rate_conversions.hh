// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// SubsTK : phylogenetic analysis and simulation library
//
//   http://www.phrap.org
//
//
// Copyright (C) 2007 Christopher T Saunders (ctsa@u.washington.edu)
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

// $Id: bioseq_util_flux_rate_conversions.hh 1085 2008-01-11 20:21:55Z ctsa $

/// \file 

#include "bioseq_util_pdf_conversions.h"
#include "../math/matrix_util.h"
#include "../math/prob_util.h"




// get w->s/s->w ratio from nuc rates and bg distro
//
template <typename FloatType1,
          typename FloatType2>
FloatType1
get_ws_ratio(const FloatType1 rates[NUC::SIZE*NUC::SIZE],
             const FloatType2 distro[]){

  FloatType1 flux[NUC::SIZE*NUC::SIZE];
  std::copy(rates,rates+NUC::SIZE*NUC::SIZE,flux);
  rates_to_flux_inplace(NUC::SIZE,distro,flux);

  const FloatType1 pws_numer =
    flux[NUC::C + NUC::A*NUC::SIZE]+
    flux[NUC::C + NUC::T*NUC::SIZE]+
    flux[NUC::G + NUC::A*NUC::SIZE]+
    flux[NUC::G + NUC::T*NUC::SIZE];
  const FloatType1 pws_denom = distro[NUC::A]+distro[NUC::T];

  const FloatType1 pws = pws_numer/pws_denom;

  const FloatType1 psw_numer =
    flux[NUC::A + NUC::C*NUC::SIZE]+
    flux[NUC::A + NUC::G*NUC::SIZE]+
    flux[NUC::T + NUC::C*NUC::SIZE]+
    flux[NUC::T + NUC::G*NUC::SIZE];
  const FloatType1 psw_denom = distro[NUC::C]+distro[NUC::G];

  const FloatType1 psw = psw_numer/psw_denom;

  return pws/psw;
}




// Get kappa terms from nuc rates and bg distro
//
template <typename FloatType1,
          typename FloatType2>
void
get_si_sv(FloatType1& transi,
          FloatType1& transv,
          const FloatType1 rates[NUC::SIZE*NUC::SIZE],
          const FloatType2 distro[]){

  // divide substitutions into prob of transition/transversion/conservation
  // ...make sure these numbers add up to one
  //
  // note that denominator for all of these probs is P(A)+P(C)+P(G)+P(T)=1.
  //

  FloatType1 flux[NUC::SIZE*NUC::SIZE];
  std::copy(rates,rates+NUC::SIZE*NUC::SIZE,flux);
  rates_to_flux_inplace(NUC::SIZE,distro,flux);

  transi =
    flux[NUC::G + NUC::A*NUC::SIZE]+
    flux[NUC::A + NUC::G*NUC::SIZE]+
    flux[NUC::T + NUC::C*NUC::SIZE]+
    flux[NUC::C + NUC::T*NUC::SIZE];

  transv =
    flux[NUC::C + NUC::A*NUC::SIZE]+
    flux[NUC::T + NUC::A*NUC::SIZE]+
    flux[NUC::A + NUC::C*NUC::SIZE]+
    flux[NUC::G + NUC::C*NUC::SIZE]+
    flux[NUC::C + NUC::G*NUC::SIZE]+
    flux[NUC::T + NUC::G*NUC::SIZE]+
    flux[NUC::A + NUC::T*NUC::SIZE]+
    flux[NUC::G + NUC::T*NUC::SIZE];

  transv /= 2.;
}




template <typename FloatType1,
          typename FloatType2>
void
rates_to_flux_inplace(const unsigned state_size,
                      const FloatType1* bg_pdistro,
                      FloatType2* rf_matrix){

  for(unsigned s1(0);s1<state_size;++s1){
    for(unsigned s2(0);s2<state_size;++s2){
      rf_matrix[s2+s1*state_size] *= bg_pdistro[s1];
    }
  }
}




template <typename FloatType1,
          typename FloatType2>
void
flux_to_rates_inplace(const unsigned state_size,
                      const FloatType1* bg_pdistro,
                      FloatType2* rf_matrix){

  for(unsigned s1(0);s1<state_size;++s1){
    const FloatType1 scale(1./bg_pdistro[s1]);
    for(unsigned s2(0);s2<state_size;++s2){
      rf_matrix[s2+s1*state_size] *= scale;
    }
  }
}




template <typename FloatType>
void
site_model_flux_to_nscodon_flux(const SITE_MODEL::index_t sm,
                                const FloatType* flux,
                                FloatType* nscodon_flux){

  const unsigned state_size(SITE_MODEL::state_size(sm));

  array_zero(nscodon_flux,NSCODON::SIZE*NSCODON::SIZE);

  for(unsigned s1(0);s1<state_size;++s1){
    const NSCODON::index_t s1nc(SITE_MODEL::decode_nscodon(sm,s1));
    for(unsigned s2(s1+1);s2<state_size;++s2){
      const NSCODON::index_t s2nc(SITE_MODEL::decode_nscodon(sm,s2));
      if(s1nc == s2nc) continue;

      nscodon_flux[s1nc+s2nc*NSCODON::SIZE] += flux[s1+s2*state_size];
      nscodon_flux[s2nc+s1nc*NSCODON::SIZE] += flux[s2+s1*state_size];
    }
  }
  matrix_balance_diagonal(nscodon_flux,NSCODON::SIZE);
}




// static
template <typename FloatType>
void
nscodon_flux_to_nsaa_flux(const FloatType nscodon_flux[NSCODON::SIZE*NSCODON::SIZE],
                          FloatType nsaa_flux[NSAA::SIZE*NSAA::SIZE]){

  array_zero(nsaa_flux,NSAA::SIZE*NSAA::SIZE);

  for(unsigned c1(0);c1<NSCODON::SIZE;++c1){
    const NSAA::index_t c1aa(codon_trans_known(static_cast<NSCODON::index_t>(c1)));
    for(unsigned c2(c1+1);c2<NSCODON::SIZE;++c2){
      const NSAA::index_t c2aa(codon_trans_known(static_cast<NSCODON::index_t>(c2)));
      if(c1aa == c2aa) continue;

      nsaa_flux[c1aa+c2aa*NSAA::SIZE] +=  nscodon_flux[c1+c2*NSCODON::SIZE];
      nsaa_flux[c2aa+c1aa*NSAA::SIZE] +=  nscodon_flux[c2+c1*NSCODON::SIZE];
    }
  }
  matrix_balance_diagonal(nsaa_flux,NSAA::SIZE);
}




/// \brief nsaa_flux is aa exchange rate (w/ selection) * P(from state)
///
template <typename FloatType>
void
neutral_nscodon_flux_to_nsaa_flux(const FloatType nscodon_flux[NSCODON::SIZE*NSCODON::SIZE],
                                  const FloatType nsaa_selection[NSAA::SIZE*NSAA::SIZE],
                                  FloatType nsaa_flux[NSAA::SIZE*NSAA::SIZE]){

  nscodon_flux_to_nsaa_flux(nscodon_flux,nsaa_flux);
  for(unsigned i(0);i<NSAA::SIZE*NSAA::SIZE;++i){ nsaa_flux[i] *= nsaa_selection[i]; }
  matrix_balance_diagonal(nsaa_flux,NSAA::SIZE);
}



#if 0
/// \brief nsaa_rate is aa exchange rate (w/ selection)
///
template <typename FloatType>
void
nscodon_flux_to_nsaa_rate(const FloatType nscodon_flux[NSAA::SIZE*NSAA::SIZE],
                          const FloatType nscodon_distro[NSCODON::SIZE]
FloatType nsaa_rate[NSAA::SIZE*NSAA::SIZE],

                                   const FloatType nsaa_selection[NSAA::SIZE*NSAA::SIZE],
                                   const FloatType nscodon_distro[NSCODON::SIZE]){

  get_nsaa_flux_from_neutral_nscodon(nsaa_rate,nscodon_rates,nsaa_selection,nscodon_distro);

  FloatType nsaa_distro[NSAA::SIZE];
  nscodon_pdf_2_nsaa_pdf(nscodon_distro,nsaa_distro);

  for(unsigned aa1(0); aa1<NSAA::SIZE; ++aa1){
    for(unsigned aa2(0); aa2<NSAA::SIZE; ++aa2){
      nsaa_rate[aa2+aa1*NSAA::SIZE] /= nsaa_distro[aa1];
    }
  }
}
#endif



// static
template <typename FloatType>
FloatType
neutral_nsaa_flux_to_dn_ds(const FloatType neutral_nsaa_flux[NSAA::SIZE*NSAA::SIZE],
                           const FloatType nsaa_selection[NSAA::SIZE*NSAA::SIZE]){

  // neutral nonsyn/syn expectation:
  FloatType lca_kn(0.);
  for(unsigned i(0);i<NSAA::SIZE;++i){
    for(unsigned j(0);j<NSAA::SIZE;++j){
      if(i==j) continue;
      lca_kn += neutral_nsaa_flux[j+i*NSAA::SIZE];
    }
  }

  // selected nonsyn/syn expectation:
  FloatType select_lca_kn(0.);
  for(unsigned i(0);i<NSAA::SIZE;++i){
    for(unsigned j(0);j<NSAA::SIZE;++j){
      if(i==j) continue;
      select_lca_kn += neutral_nsaa_flux[j+i*NSAA::SIZE]*nsaa_selection[j+i*NSAA::SIZE];
    }
  }

  return select_lca_kn/lca_kn;
}




template <typename FloatType>
FloatType
neutral_nscodon_flux_to_dn_ds(const FloatType neutral_nscodon_flux[NUC::SIZE*NUC::SIZE],
                              const FloatType nsaa_selection[NSAA::SIZE*NSAA::SIZE]){

  FloatType neutral_nsaa_flux[NSAA::SIZE*NSAA::SIZE];
  nscodon_flux_to_nsaa_flux(neutral_nscodon_flux,neutral_nsaa_flux);

  return neutral_nsaa_flux_to_dn_ds(neutral_nsaa_flux,nsaa_selection);
}





// static
template <typename FloatType>
void
neutral_nsaa_flux_to_nsaa_dn_ds(const FloatType neutral_nsaa_flux[NSAA::SIZE*NSAA::SIZE],
                                const FloatType nsaa_selection[NSAA::SIZE*NSAA::SIZE],
                                FloatType to_nsaa_dn_ds[NSAA::SIZE],
                                FloatType from_nsaa_dn_ds[NSAA::SIZE]){

  for(unsigned i(0);i<NSAA::SIZE;++i){
    // neutral expected flux out of amino acid i:
    FloatType from_kn(0.);
    for(unsigned j(0);j<NSAA::SIZE;++j){
      if(i==j) continue;
      from_kn += neutral_nsaa_flux[j+i*NSAA::SIZE];
    }

    // selected flux out of amino acid i
    FloatType select_from_kn(0.);
    for(unsigned j(0);j<NSAA::SIZE;++j){
      if(i==j) continue;
      select_from_kn += neutral_nsaa_flux[j+i*NSAA::SIZE]*nsaa_selection[j+i*NSAA::SIZE];
    }

    from_nsaa_dn_ds[i] = select_from_kn/from_kn;
  }

  for(unsigned j(0);j<NSAA::SIZE;++j){
    // neutral expected flux into amino acid j:
    FloatType to_kn(0.);
    for(unsigned i(0);i<NSAA::SIZE;++i){
      if(i==j) continue;
      to_kn += neutral_nsaa_flux[j+i*NSAA::SIZE];
    }

    // selected flux into amino acid j
    FloatType select_to_kn(0.);
    for(unsigned i(0);i<NSAA::SIZE;++i){
      if(i==j) continue;
      select_to_kn += neutral_nsaa_flux[j+i*NSAA::SIZE]*nsaa_selection[j+i*NSAA::SIZE];
    }

    to_nsaa_dn_ds[j] = select_to_kn/to_kn;
  }
}




template <typename FloatType>
void
neutral_nscodon_flux_to_nsaa_dn_ds(const FloatType neutral_nscodon_flux[NUC::SIZE*NUC::SIZE],
                                   const FloatType nsaa_selection[NSAA::SIZE*NSAA::SIZE],
                                   FloatType to_nsaa_dn_ds[NSAA::SIZE],
                                   FloatType from_nsaa_dn_ds[NSAA::SIZE]){

  FloatType neutral_nsaa_flux[NSAA::SIZE*NSAA::SIZE];
  nscodon_flux_to_nsaa_flux(neutral_nscodon_flux,neutral_nsaa_flux);

  neutral_nsaa_flux_to_nsaa_dn_ds(neutral_nsaa_flux,nsaa_selection,to_nsaa_dn_ds,from_nsaa_dn_ds);
}

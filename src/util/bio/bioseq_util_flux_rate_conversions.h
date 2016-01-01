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

// $Id: bioseq_util_flux_rate_conversions.h 1085 2008-01-11 20:21:55Z ctsa $

/// \file 
///
/// \brief a variety of flux/rate and parameter summary tools. file is probably misnamed
///



//   summary of condensing rate matrices to a smaller
//   state subset:
//
//   for A,B in S
//   and X,Y in T
//
//   get P(S|T) given P(A|X),P(A|Y),P(B|X),P(B|Y) :
//
//   P(to=S|from=T) = P(to=S,from=T)/P(from=T)
//
//   P(to=S,from=T) = P(A|X)P(X)+P(A|Y)P(Y)+P(B|X)P(X)+P(B|Y)P(Y)
//   P(from=T) = P(X)+P(Y)
//


#ifndef __BIOSEQ_UTIL_FLUX_RATE_CONVERSIONS_H
#define __BIOSEQ_UTIL_FLUX_RATE_CONVERSIONS_H

#include "bioseq_util.h"

/// \brief Get the w->s/s->w ratio from rates and ancestral distro
///
template <typename FloatType1,
          typename FloatType2>
FloatType1
get_ws_ratio(const FloatType1 rates[NUC::SIZE*NUC::SIZE],
             const FloatType2 distro[]);

/// \brief Get transition/transversion/conservation rates
///
template <typename FloatType1,
          typename FloatType2>
void
get_si_sv(FloatType1& transi,
          FloatType1& transv,
          const FloatType1 rates[NUC::SIZE*NUC::SIZE],
          const FloatType2 distro[]);

template <typename FloatType1,
          typename FloatType2>
void
rates_to_flux_inplace(const unsigned state_size,
                      const FloatType1* bg_pdistro,
                      FloatType2* rf_matrix);

template <typename FloatType1,
          typename FloatType2>
void
flux_to_rates_inplace(const unsigned state_size,
                      const FloatType1* bg_pdistro,
                      FloatType2* rf_matrix);

template <typename FloatType>
void
site_model_flux_to_nscodon_flux(const SITE_MODEL::index_t sm,
                                const FloatType* flux,
                                FloatType* nscodon_flux);

/// \brief nsaa_flux is aa exchange rate (w/ selection) * P(from state)
///
template <typename FloatType>
void
neutral_nscodon_flux_to_nsaa_flux(const FloatType nscodon_flux[NSCODON::SIZE*NSCODON::SIZE],
                                  const FloatType nsaa_selection[NSAA::SIZE*NSAA::SIZE],
                                  FloatType nsaa_flux[NSAA::SIZE*NSAA::SIZE]);

/// \brief condense nsaa dN/dS matrix to global dN/dS:
template <typename FloatType>
FloatType
neutral_nscodon_flux_to_dn_ds(const FloatType neutral_nscodon_flux[NUC::SIZE*NUC::SIZE],
                              const FloatType nsaa_selection[NSAA::SIZE*NSAA::SIZE]);

/// \brief condense nsaa dN/dS matrix to aa to/from dN/dS:
template <typename FloatType>
void
neutral_nscodon_flux_to_nsaa_dn_ds(const FloatType neutral_nscodon_flux[NUC::SIZE*NUC::SIZE],
                                   const FloatType nsaa_selection[NSAA::SIZE*NSAA::SIZE],
                                   FloatType to_nsaa_dn_ds[NSAA::SIZE],
                                   FloatType from_nsaa_dn_ds[NSAA::SIZE]);

#include "bioseq_util_flux_rate_conversions.hh"

#endif

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

// $Id: rate_gtor_nscodon_base_util.h 1145 2008-02-28 20:07:28Z ctsa $

/// \file

#ifndef __RATE_GTOR_NSCODON_BASE_UTIL_H
#define __RATE_GTOR_NSCODON_BASE_UTIL_H

#include "rate_gtor_nscodon_base.h"
#include "subs_ml_types.h"
#include "util/bio/bioseq_util.h"

template <typename FloatType>
void
nscodon_flux_no_aa_selection(const rate_gtor_nscodon_base& r,
                             const unsigned cat_no,
                             const unsigned branch_cat_set,
                             FloatType* nscodon_flux);

/// \brief expected selection on mutation given partial codon context
///
/// given a nucleotide mutation at a single position in a codon,
/// calculate the average selection against this mutations for all
/// possible codons formed by the remaining nucleotides
///
smlfloat
average_codon_selection(const rate_gtor_nscodon_base& r,
                        const NUC::index_t* nuc,
                        const NUC::index_t mut_nuc,
                        const unsigned mut_pos,
                        const rates_func_options& opt,
                        const prob_t* nscodon_pdistro);

/// \brief produce rate matrix for site-model sm from r
///
template <typename FloatType>
void
rates_nscodon_context(const rate_gtor_nscodon_base& r,
                      const SITE_MODEL::index_t sm,
                      FloatType rates[],
                      const rates_func_options& opt);

/// \brief produce rate matrix for native site-model from r
///
template <typename FloatType>
void
rates_nscodon_context(const rate_gtor_nscodon_base& r,
                      FloatType rates[],
                      const rates_func_options& opt){

  rates_nscodon_context(r,r.site_model(),rates,opt);
}


#include "rate_gtor_nscodon_base_util.hh"

#endif

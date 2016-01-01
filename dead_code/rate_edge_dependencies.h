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

// $Id: rate_edge_dependencies.h 1099 2008-01-21 22:54:47Z ctsa $
//

/// \file

#ifndef __RATE_EDGE_DEPENDENCIES_H
#define __RATE_EDGE_DEPENDENCIES_H

#include "rate_gtor_options.h"
#include "subs_ml_types.h"
#include "util/bio/bioseq_util.h"

#if 0

void
shape_dinuc_nscodon(const prob_t* x_nscodon,
                    prob_t* x_dinuc);

#endif

#ifdef BRANCH_SPECIFIC_OBS
struct site_model_edge_pdistros {

  void ropt_hookup(rates_func_options& ropt) const;

  prob_t nuc_pdistro_overedge_5p[NUC::SIZE];
  prob_t nuc_pdistro_overedge_3p[NUC::SIZE];
  prob_t nuc_pdistro_overedge_5p_cond_on_edge[NUC::SIZE*NUC::SIZE];
  prob_t nuc_pdistro_overedge_3p_cond_on_edge[NUC::SIZE*NUC::SIZE];
  prob_t nscodon_pdistro[NSCODON::SIZE];
};

void
set_site_model_edge_pdistros(const SITE_MODEL::index_t sm,
                             const prob_t* sm_pdistro,
                             site_model_edge_pdistros& smp);
#endif



// old site edge handling system: used the stationary distro
// conditioned on the (known) edge nuc to get the distribution of the
// unknown edge nuc -- it works pretty well when the root is locked to
// the present day state distro
//

#if 0
/// solve for dependency between an "edge" nucleotide, at nucleotide
/// position=edge_pos in the independent site, and another "overedge"
/// nucleotide either before (overedge==-1) or after (overedge==+1) the
/// edge nucleotide.
///
void
get_rate_edge_dependencies_noncoding_unidirectional_context(const rate_gtor_nuc_base& r,
                                                            const unsigned site_category,
                                                            const unsigned group_category,
                                                            const int overedge_dir,
                                                            prob_t* bg_pdf_nuc_overedge_cond_on_edge,
                                                            prob_t* bg_pdf_nuc_edge_cond_on_overedge = 0);

/// solve for dependency between an "edge" nucleotide, and another
/// "overedge" nucleotide either before (overedge_dir==-1) or after
/// (overedge_dir==+1) the edge nucleotide. When context is
/// bidirectional, this is done by stepping thru all values for the
/// "underedge" (edge-overedge_dir) nuc, and building a rate matrix for
/// the exchange of the edge and overedge nuc conditioned on the
/// underedge value. The stationary distro of this rate matrix is
/// solved, yielding the edge,overedge joint, which is then transformed
/// to the appropriate conditional distribution
///
void
get_rate_edge_dependencies_noncoding_bidirectional_context(const rate_gtor_nuc_base& r,
                                                           const unsigned site_category,
                                                           const unsigned group_category,
                                                           const int overedge_dir,
                                                           prob_t bg_pdf_nuc_overedge_cond_on_edge_underedge[NUC::SIZE][NUC::SIZE][NUC::SIZE]);

#ifdef STAT_EDGE_CODING
/// solve for dependency between an "edge" nucleotide in the codon,
/// (edge_pos), and another "overedge" nucleotide either before
/// (overedge_dir==-1) or after (overedge_dir==+1) it.
///
///  throws subs_ml_error
void
get_rate_edge_dependencies_coding(const rate_gtor_nscodon_base& r,
                                  const rates_func_options& opt,
                                  const unsigned edge_pos,
                                  const int overedge_dir,
                                  prob_t* bg_pdf_nuc_overedge_cond_on_edge,
                                  prob_t* bg_pdf_nuc_edge_cond_on_overedge = 0);
#endif
#endif



#endif

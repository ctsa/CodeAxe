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

// $Id: rates_func_options.h 1107 2008-01-25 02:40:20Z ctsa $

/// \file

#ifndef __RATES_FUNC_OPTIONS_H
#define __RATES_FUNC_OPTIONS_H

#include "subs_ml_types.h"



struct rates_func_options_base {
  explicit
  rates_func_options_base(unsigned c,
                          bool saa = false,
                          bool snonaa = false,
                          bool is_use_sm = false,
                          unsigned sm = 0)
    : cat(c),
      is_skip_aa_selection(saa),is_skip_nonaa_selection(snonaa),
      is_use_submodel(is_use_sm), submodel_no(sm),
      nuc_pdistro_overedge_5p(0),
      nuc_pdistro_overedge_3p(0),
      nuc_pdistro_overedge_5p_cond_on_edge(0),
      nuc_pdistro_overedge_3p_cond_on_edge(0),
      nscodon_pdistro(0) {}

  unsigned cat;
  bool is_skip_aa_selection;
  bool is_skip_nonaa_selection;
  bool is_use_submodel;
  int submodel_no;
  const prob_t* nuc_pdistro_overedge_5p;// [NUC::SIZE]
  const prob_t* nuc_pdistro_overedge_3p;// [NUC::SIZE]
  const prob_t* nuc_pdistro_overedge_5p_cond_on_edge;// [NUC::SIZE*NUC::SIZE]
  const prob_t* nuc_pdistro_overedge_3p_cond_on_edge;// [NUC::SIZE*NUC::SIZE]
  const prob_t* nscodon_pdistro;// [NSCODON::SIZE]
};



struct rates_func_options : public rates_func_options_base {

  typedef rates_func_options_base base_t;

  rates_func_options(const rates_func_options_base& rfob,
                     const unsigned bcs)
    :  base_t(rfob), branch_cat_set(bcs) {}

  unsigned branch_cat_set;
};


#endif

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

// $Id: rate_gtor_dinuc.h 1107 2008-01-25 02:40:20Z ctsa $

/// \file

#ifndef __RATE_GTOR_DINUC_H
#define __RATE_GTOR_DINUC_H


#include "condition_func_dinuc.h"
#include "rate_gtor_nuc_base.h"
#include "rate_gtor_nuc_base_util.h"

/// \brief double nuc site conditioned to single nuc
///
struct rate_gtor_dinuc : public rate_gtor_nuc_base {

  typedef rate_gtor_nuc_base base_t;
  typedef rate_gtor_dinuc self_t;

  rate_gtor_dinuc(const rate_gtor_options& ropt,
                  const rate_gtor_nuc_options& nopt,
                  const rate_gtor_sml_share& smls)
    : base_t(ropt,nopt,smls) {}

  rate_gtor_dinuc(const self_t& s,
                  const rate_gtor_sml_share& smls) : base_t(s,smls) {}

  virtual
  rate_gtor* clone(const rate_gtor_sml_share& smls) const { return new self_t(*this,smls); }

  virtual
  unsigned state_size_conditioned() const { return NUC::SIZE; }

  virtual
  condition_func* condition_func_factory() const { return new condition_func_dinuc;}

  ///////////////////////////////////////////////////////////////////////////
private:

  virtual
  void state_pdistro_to_nuc_pdistro(const prob_t* state_pdistro,
                                    prob_t* nuc_pdistro) const;
};


#endif

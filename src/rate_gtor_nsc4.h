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

// $Id: rate_gtor_nsc4.h 1107 2008-01-25 02:40:20Z ctsa $

/// \file

#ifndef __RATE_GTOR_NSC4_H
#define __RATE_GTOR_NSC4_H


#ifdef USE_BIDIR_NSC4
#include "condition_func_nsc4_bidir.h"
#else
#include "condition_func_nsc4.h"
#endif
#include "rate_gtor_nscodon_base.h"
#include "util/general/die.h"

/// \brief rate generator for C4+/C4- states
///
struct rate_gtor_nsc4 : public rate_gtor_nscodon_base {

  typedef rate_gtor_nscodon_base base_t;
  typedef rate_gtor_nsc4 self_t;

  rate_gtor_nsc4(const rate_gtor_options& ropt,
                 const rate_gtor_nuc_options& nopt,
                 const rate_gtor_nscodon_options& copt,
                 const rate_gtor_sml_share& smls)
    : base_t(ropt,nopt,copt,smls) {}

  rate_gtor_nsc4(const self_t& s,
                 const rate_gtor_sml_share& smls) : base_t(s,smls) {}

  ////////////////////////////////////////////////////////////////////////////
private:

  virtual
  rate_gtor* clone(const rate_gtor_sml_share& smls) const { return new self_t(*this,smls); }

  virtual
  unsigned state_size_conditioned() const { return NSCODON::SIZE; }

  virtual
  condition_func* condition_func_factory() const {
#ifdef USE_BIDIR_NSC4
    return new condition_func_nsc4_bidir(site_model());
#else
    return new condition_func_nsc4();
#endif
  }
};


#endif

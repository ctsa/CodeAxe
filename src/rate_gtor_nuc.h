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

// $Id: rate_gtor_nuc.h 1107 2008-01-25 02:40:20Z ctsa $

/// \file

#ifndef __RATE_GTOR_NUC_H
#define __RATE_GTOR_NUC_H

#include "rate_gtor_nuc_base.h"
#include "util/bio/bioseq_util.h"

#include <iosfwd>


/// \brief rate generator for nuc model
///
struct rate_gtor_nuc : public rate_gtor_nuc_base {

  typedef rate_gtor_nuc_base base_t;
  typedef rate_gtor_nuc self_t;

  rate_gtor_nuc(const rate_gtor_options& ropt,
                const rate_gtor_nuc_options& nopt,
                const rate_gtor_sml_share& smls)
    : base_t(ropt,nopt,smls) {}

  rate_gtor_nuc(const self_t& s,
                const rate_gtor_sml_share& smls) : base_t(s,smls) {}

  virtual
  rate_gtor* clone(const rate_gtor_sml_share& smls) const { return new self_t(*this,smls);}

private:
  virtual
  void state_pdistro_to_nuc_pdistro(const prob_t* state_pdistro,
                                    prob_t* nuc_pdistro) const {
    for(unsigned i(0);i<NUC::SIZE;++i) { nuc_pdistro[i] = state_pdistro[i]; }
  }
};

#endif

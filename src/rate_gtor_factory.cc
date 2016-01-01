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

// $Id: rate_gtor_factory.cc 1107 2008-01-25 02:40:20Z ctsa $

/// \file

#include "rate_gtor_binary.h"
#include "rate_gtor_dinuc.h"
#include "rate_gtor_factory.h"
#include "rate_gtor_nsc4.h"
#include "rate_gtor_nsc5.h"
#include "rate_gtor_nscodon.h"
#include "rate_gtor_nuc.h"
#include "rate_gtor_trinuc.h"
#include "util/general/log.h"



rate_gtor*
rate_gtor_factory(const rate_gtor_sml_share& smls,
                  const rate_gtor_options& ropt,
                  const rate_gtor_nuc_options& nopt,
                  const rate_gtor_nscodon_options& copt,
                  const C5_APPROX::index_t c5t){

  using namespace RATE_GTOR_MODEL;

  const index_t rgm(smls.rate_gtor_model());

  switch(rgm){
  case BINARY:                  return new rate_gtor_binary(ropt,smls);
  case RATE_GTOR_MODEL::NUC:    return new rate_gtor_nuc(ropt,nopt,smls);
  case RATE_GTOR_MODEL::CODON:  return new rate_gtor_nscodon(ropt,nopt,copt,smls);
  case C4PRE:
  case C4POST:                  return new rate_gtor_nsc4(ropt,nopt,copt,smls);
  case C5:
  case C5PRE:
  case C5POST:                  return new rate_gtor_nsc5(ropt,nopt,copt,c5t,smls);
  case RATE_GTOR_MODEL::DINUC:  return new rate_gtor_dinuc(ropt,nopt,smls);
  case RATE_GTOR_MODEL::TRINUC: return new rate_gtor_trinuc(ropt,nopt,smls);
  case NONE:                    return static_cast<rate_gtor*>(0);
  default:
    log_os << "ERROR:: unknown model type in call to rate_gtor_factory: "
           << static_cast<int>(rgm) << "\n";
    abort();
    return static_cast<rate_gtor*>(0);
  }
}



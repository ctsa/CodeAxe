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

// $Id: prob_util_io.h 743 2007-08-14 15:47:12Z ctsa $

/// \file
///
/// \brief i/o functions for small discrete prob distros
///

#ifndef __PROB_UTIL_IO_H
#define __PROB_UTIL_IO_H

#include <iosfwd>


/// \brief simple debugging output
///
template <typename FloatType>
void
pdistro_dump(const FloatType pn[],
             const unsigned N,
             std::ostream& os,
             const char* label=0);

/// \brief output w/ syms:
///
template <typename FloatType,
          typename SymType>
void
pdistro_report(const FloatType* pdf,
               const unsigned N,
               const SymType* syms,
               std::ostream& os);

#if 0
/// \brief pretty print a transition (stochastic) matrix
template <typename FloatType>
void
report_transition_prob(const FloatType* tprob, // dim*dim
                       const char* label,
                       const char* syms,
                       std::ostream& os,
                       const unsigned dim);


template <typename FloatType>
void
report_transition_prob_flux(const FloatType* tprob,       // dim*dim
                            const FloatType* from_distro, // dim
                            const char* label,
                            const char* syms,
                            std::ostream& os,
                            const unsigned dim);
#endif


#include "prob_util_io.hh"

#endif

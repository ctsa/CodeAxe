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

// $Id: bioseq_util_io.h 585 2007-06-12 20:16:09Z ctsa $

/// \file 

#ifndef __BIOSEQ_UTIL_IO_H
#define __BIOSEQ_UTIL_IO_H


#include "bioseq_util.h"

#include <ostream>


template <typename FloatType>
void
print_site_model_distro(const FloatType pdf[],
                        const SITE_MODEL::index_t sm,
                        std::ostream& os){

  os << SITE_MODEL::label(sm) << " distro:\n";

  const unsigned ss(SITE_MODEL::state_size(sm));
  for(unsigned i(0);i<ss;++i){
    os << SITE_MODEL::print(sm,i) << ": " << pdf[i] << "\n";
  }
  os << "\n";
}


#endif

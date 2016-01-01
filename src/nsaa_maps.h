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

// $Id: nsaa_maps.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __NSAA_MAPS_H
#define __NSAA_MAPS_H

#include "util/bio/bioseq_util.h"

namespace NSAA_MAPS {

  namespace LEHNINGER {
    enum index_t {
      ALIPHATIC,
      AROMATIC,
      POLAR,
      POSITIVE,
      NEGATIVE,
      SIZE
    };

    extern const char syms[SIZE];

    const unsigned*
    get_reduction_map();
  }

  namespace LANGE {
    enum index_t {
      NONPOLAR,
      POSITIVE,
      NEGATIVE,
      CYS,
      SIZE
    };

    extern const char syms[SIZE];

    const unsigned*
    get_reduction_map();
  }

  namespace HP {
    enum index_t {
      H,
      P,
      SIZE
    };

    extern const char syms[SIZE];

    const unsigned*
      get_reduction_map();
  }

  namespace HPO {
    enum index_t {
      H,
      P,
      O,
      SIZE
    };

    extern const char syms[SIZE];

    const unsigned*
      get_reduction_map();
  }

  namespace BETA_BRANCH {
    enum index_t {
      OTHER,
      BETA_BRANCH,
      SIZE
    };

    extern const char syms[SIZE];

    const unsigned*
    get_reduction_map();
  }

  namespace HELIX_FORM {
    enum index_t {
      OTHER,
      HELIX_FORM,
      SIZE
    };

    extern const char syms[SIZE];

    const unsigned*
    get_reduction_map();
  }

  namespace AROMATIC {
    enum index_t {
      OTHER,
      AROMATIC,
      SIZE
    };

    extern const char syms[SIZE];

    const unsigned*
    get_reduction_map();
  }

  namespace CYS {
    enum index_t {
      OTHER,
      CYS,
      SIZE
    };

    extern const char syms[SIZE];

    const unsigned*
    get_reduction_map();
  }

  namespace GLY {
    enum index_t {
      OTHER,
      GLY,
      SIZE
    };

    extern const char syms[SIZE];

    const unsigned*
    get_reduction_map();
  }


  namespace PRO {
    enum index_t {
      OTHER,
      PRO,
      SIZE
    };

    extern const char syms[SIZE];

    const unsigned*
    get_reduction_map();
  }


  struct converter {

    static
    HP::index_t
    aa_to_hp(const NSAA::index_t a){
      return static_cast<HP::index_t>(hpmap[a]);
    }

    static const unsigned* hpmap;
  };

}


#endif

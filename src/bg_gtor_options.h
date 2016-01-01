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

// $Id: bg_gtor_options.h 1117 2008-01-28 19:48:39Z ctsa $

/// \file

#ifndef __BG_GTOR_OPTIONS_H
#define __BG_GTOR_OPTIONS_H



namespace SUBS_RATE_BG_MODEL {
  enum index_t { OBS_AVG,
                 FULL,
                 CODON_DINUC
  };
}



struct bg_gtor_options {

  bg_gtor_options() :
    bg(SUBS_RATE_BG_MODEL::OBS_AVG) {}

  SUBS_RATE_BG_MODEL::index_t bg;
};



#endif

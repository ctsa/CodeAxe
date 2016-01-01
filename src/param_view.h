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

// $Id: param_view.h 744 2007-08-14 18:09:31Z ctsa $

/// \file

#ifndef __PARAM_VIEW_H
#define __PARAM_VIEW_H


namespace PARAM_VIEW {

  enum index_t {
    ALL,
    TRAINABLE,     // trainable only, untransformed
    MIN_PARAM,     // minimization transformed trainable + pseudo parameters
    INDY_MIN_PARAM // minimization transformed true independent parameters derived from trainable (for gradient minimizers)
  };
}


#endif

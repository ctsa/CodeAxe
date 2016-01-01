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

// $Id: minfunc_interface.h 1033 2007-12-03 22:25:57Z ctsa $

/// \file

#ifndef __MINFUNC_INTERFACE_H
#define __MINFUNC_INTERFACE_H

/// \brief min function interface for line minimizers
///
template <typename FloatType>
struct minfunc_1d_interface {
  virtual ~minfunc_1d_interface() {}

  virtual FloatType val(const FloatType x) const  = 0;

  virtual bool is_val_computable(const FloatType) const { return true; }
};


/// \brief min function interface for multidimensional minimizers
///
template <typename FloatType>
struct minfunc_interface {
  virtual ~minfunc_interface() {}

  virtual unsigned dim() const = 0;

  virtual FloatType val(const FloatType* v) = 0;

  // this is not guaranteed to be called before val(), it only gives
  // minimizers the option to handle out-of-bounds parameter space
  // more gracefully than would be possible with a penalized score
  //
  virtual bool is_val_computable(const FloatType*) { return true; }
};


template <typename FloatType>
struct minfunc_gradient_interface : public minfunc_interface<FloatType> {

  virtual
  FloatType dval(const FloatType* v,
                 FloatType* dv) = 0;
};


#endif

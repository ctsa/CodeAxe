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

// $Id: deriv_util.h 1065 2007-12-14 21:49:28Z ctsa $

/// \file

#ifndef __DERIVE_UTIL_H
#define __DERIVE_UTIL_H

#include <cmath>



template <typename FloatType,
          typename Func>
FloatType
forward_deriv_estimate(const FloatType x,
                       const FloatType fx,
                       Func& f,
                       const FloatType delta = 1e-5){

  const FloatType deltamin(delta*1e-2);
  const FloatType xdelta(std::max(std::fabs(x*delta),deltamin));

  const FloatType x_hi(x+xdelta);
  const FloatType fx_hi(f(x_hi));

  return (fx_hi-fx)/(xdelta);
}



template <typename FloatType,
          typename Func>
FloatType
centered_deriv_estimate(const FloatType x,
                        Func& f,
                        const FloatType delta = 1e-5){

  const FloatType deltamin(delta*1e-2);
  const FloatType xdelta(std::max(std::fabs(x*delta),deltamin));

  const FloatType xhi(x+xdelta);
  const FloatType xlo(x-xdelta);

  const FloatType fxhi(f(xhi));
  const FloatType fxlo(f(xlo));

  return (fxhi-fxlo)/(2*xdelta);
}




// (df^2)/(dx1,dx2)
//
template <typename FloatType,
          typename Func>
FloatType
centered_2partial_deriv_estimate(const FloatType x1,
                                 const FloatType x2,
                                 Func& f,
                                 const FloatType delta = 1e-5,
                                 const FloatType deltamin_factor = 1e-2){

  const FloatType deltamin(delta*deltamin_factor);
  const FloatType x1delta(std::max(std::fabs(x1*delta),deltamin));
  const FloatType x2delta(std::max(std::fabs(x2*delta),deltamin));

  const FloatType x1hi(x1+x1delta);
  const FloatType x1lo(x1-x1delta);
  const FloatType x2hi(x2+x2delta);
  const FloatType x2lo(x2-x2delta);

  const FloatType fx1hi_x2hi(f(x1hi,x2hi));
  const FloatType fx1lo_x2hi(f(x1lo,x2hi));
  const FloatType fx1hi_x2lo(f(x1hi,x2lo));
  const FloatType fx1lo_x2lo(f(x1lo,x2lo));

  return (fx1hi_x2hi-fx1lo_x2hi-fx1hi_x2lo+fx1lo_x2lo)/(4.*x1delta*x2delta);
}



template <typename FloatType,
          typename Func>
FloatType
centered_2nd_deriv_estimate(const FloatType x,
                            const FloatType fx,
                            Func& f,
                            const FloatType delta = 1e-5,
                            const FloatType deltamin_factor = 1e-2){

  const FloatType deltamin(delta*deltamin_factor);
  const FloatType xdelta(std::max(std::fabs(x*delta),deltamin));

  const FloatType xhi(x+xdelta);
  const FloatType xlo(x-xdelta);

  const FloatType fxhi(f(xhi));
  const FloatType fxlo(f(xlo));

  return (fxhi-2*fx+fxlo)/(xdelta*xdelta);
}


#endif

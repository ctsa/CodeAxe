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

// $Id: minimize_conj_gradient.hh 1055 2007-12-08 04:27:31Z ctsa $

/// \file

#include "array_util.h"
#include "minimize_line.h"

#include "../general/log.h"
#include "../general/workspace.h"

#include <exception>



template <typename FloatType,
          typename MinFunc>
void
minimize_conj_gradient(FloatType* vec_current,
                       MinFunc& mf,
                       const FloatType tol,
                       FloatType& f_min,
                       unsigned& iter,
                       FloatType& final_iter_delta_f,
                       const unsigned max_iter){

  const unsigned n(mf.dim());

  workspace<FloatType> ws(4*n);

  FloatType* vec_deriv(ws.ptr());
  FloatType* vec_grad(ws.ptr()+n);
  FloatType* vec_conj_grad(ws.ptr()+2*n);
  FloatType* linespace(ws.ptr()+3*n);

  FloatType gdot(0);

  std::fill(vec_conj_grad,vec_conj_grad+n,FloatType(0.));

  // more intuitive internal name:
  FloatType& iter_delta_f(final_iter_delta_f);

  iter=0;
  while(iter<max_iter){
    iter++;

    f_min = mf.dval(vec_current,vec_deriv);

    FloatType gnum(0);
    if(iter>1){
      for(unsigned i(0);i<n;++i){
        gnum += (vec_deriv[i]+vec_grad[i])*vec_deriv[i];
      }
      gnum /= gdot;
    }

    for(unsigned i(0);i<n;++i){
      vec_grad[i] = -vec_deriv[i];
      vec_conj_grad[i] = vec_grad[i]+gnum*vec_conj_grad[i];
    }

    const FloatType f_min_prev(f_min);

    try {
      minimize_line(vec_current,vec_conj_grad,f_min,mf,f_min,linespace);
    } catch (std::exception& e) {
      log_os << "EXCEPTION: " << e.what() << "\n";
      log_os << "...caught in " << __FILE__ << ":" << __LINE__ << " gnum,gdot " << gnum << " " << gdot << "\n";
      for(unsigned i(0);i<n;++i){
        log_os << "i, x, -grad(x), conj_grad: " << i << "\t"
               << vec_current[i] << "\t" << vec_grad[i] << "\t" << vec_conj_grad[i] << "\n";
      }
      throw;
    }

    iter_delta_f=(f_min_prev-f_min);

    if(iter_delta_f <= tol) { // converging on absolute tolerance
      if(iter_delta_f<0.) throw math_util_exception("minimize_conj_grad(): value of f increased in iteration");
      break;
    }

    gdot = array_dot(vec_grad,vec_grad,n);
    if(gdot == 0.) break;
  }
}

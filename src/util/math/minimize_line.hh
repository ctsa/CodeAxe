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

// $Id: minimize_line.hh 744 2007-08-14 18:09:31Z ctsa $

/// \file

#include "array_util.h"
#include "minfunc_interface.h"
#include "minimize_1d.h"

#include <algorithm>

#ifdef DEBUG
#include "../general/log.h"
#include <ostream>
#endif

/// \brief wraps a multidimensional function to expose a 1d interface
/// for line minimization
///
template <typename FloatType,
          typename MinFunc>
struct minfunc_1d_wrapper : public minfunc_1d_interface<FloatType> {

  minfunc_1d_wrapper(MinFunc& init_mf,
                     FloatType* workspace, // 1*mf.dim()
                     const FloatType* init_v_base,
                     const FloatType* init_v_dir,
                     const FloatType init_offset = 0)
    : n(init_mf.dim()), v_base(init_v_base), v_dir(init_v_dir),
      v_tmp(workspace), mf(init_mf), offset(init_offset) {}

  bool
  is_val_computable(const FloatType x) const {
    setup_vt(x);
    return mf.is_val_computable(v_tmp);
  }

  FloatType
  val(const FloatType x) const {
    setup_vt(x);

#ifdef DEBUG
    FloatType a;
    try {
      a = mf.val(v_tmp);
    } catch(...) {
      log_os << "UNKNOWN EXCEPTION:: minfunc_1d_wrapper, call to val(v_tmp)\n";
      for(unsigned i(0);i<n;++i){
        log_os << " i,v_tmp[i]: " << i << " " << v_tmp[i] << "\n";
      }
      throw;
    }
    return a;
#else
    return mf.val(v_tmp);
#endif
  }

private:
  void
  setup_vt(const FloatType x) const {
    std::copy(v_base,v_base+n,v_tmp);
    array_scale_plus_vector(v_dir,(x-offset),v_tmp,n);
  }

public:

  unsigned n;
  const FloatType* v_base;
  const FloatType* v_dir;
  FloatType* v_tmp;
  MinFunc& mf;
  FloatType offset;  ///<< here to prevent dealing with all d hanging around 0 when initial guess is accurate
};




template <typename FloatType,
          typename MinFunc>
void
minimize_line(FloatType* vec_base,
              FloatType* vec_dir,
              FloatType f_base,
              MinFunc& mf,
              FloatType& f_min,
              FloatType* workspace,
              FloatType x_tol,
              FloatType f_tol){

  static const FloatType x_offset(1);
  static const FloatType x1(x_offset), x2(x_offset+1);

  minfunc_1d_wrapper<FloatType,MinFunc> mf1d(mf,workspace,vec_base,vec_dir,x_offset);

  FloatType x_min;

  minimize_1d(x1,x2,f_base,mf1d,x_min,f_min,x_tol,f_tol);

  x_min -= x_offset;

  array_scale_plus_vector(vec_dir,x_min,vec_base,mf.dim());
}

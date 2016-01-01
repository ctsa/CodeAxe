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

// $Id: matrix_exp_series.hh 1055 2007-12-08 04:27:31Z ctsa $

/// \file


#include "test_float.h"

#include "../general/log.h"

#include <ostream>
#include <sstream>
#include <vector>

//#define MEXP_DEBUG

//#define ORIG_WATERMAN_METHOD


template <typename FloatType>
void
matrix_exp_series_prepdata<FloatType>::
grow(){
  if( (size+CHUNK_SIZE) > SERIES_MAX ){
    std::ostringstream oss;
    oss << "matrix_exp_series_prepdata.grow(): matrix_exp series limit reached, size: " << size;
    throw math_util_exception(oss.str().c_str());
  }

  size += CHUNK_SIZE;

  {
    const unsigned prep_size(size*dim*dim*sizeof(FloatType));
    ws1.resize_with_data(prep_size);
  }
  prep = reinterpret_cast<FloatType* const>(ws1.ptr());

  {
    const unsigned info_size(size*sizeof(iter_info));
    ws2.resize_with_data(info_size);
  }
  info = reinterpret_cast<iter_info* const>(ws2.ptr());

  for(unsigned i(size-CHUNK_SIZE);i<size;++i) info[i].is_calc = false;
}




template <typename FloatType>
void
matrix_exp_series(FloatType* out_mat,      // dim*dim
                  const FloatType* in_mat, // dim*dim
                  const FloatType in_scale,
                  const unsigned dim){

  matrix_exp_series_prepdata<FloatType> mepd(dim);

  matrix_exp_series_prep<FloatType>(mepd,in_mat);

  matrix_exp_series_scale<FloatType>(out_mat,mepd,in_scale);
}



// pre-calc the non-scaled taylor series matrices, so that different times
// for the same matrix can be quickly calculated.
//
// ...this #ifdef is just for humans:
#ifdef ORIG_WATERMAN_METHOD
// added stochastic balancing suggested by Waterman p. 369
//
// for rate matrix Q:
// 1) maxd = -min(diag(Q))
// 2) set Q' = (Q./maxd + I)
// 3) instead of exp(Qt) = sum_n(((Q.*t)^n)/(n!))
//    we calc    exp(Qt) = exp(-maxd*t) * sum_n((maxd^n)*((Q'.*t)^n)/(n!))
//
// ..this is found from Q = (Q' - I).*maxd
//                      exp(Q.*t) = exp(Q'.*maxd.*t)*exp(-I.*maxd.*t)
//
// 4) we separate out time terms so that the time independent terms can
//   be cached for multi-time calls:
//               exp(Qt) = exp(-maxd*t) * sum_n(((maxd*t)^n)*(Q'^n)/(n!))
//
// note that in matrix_exp_scale below:
// emaxd = exp(-maxd*t)
// scale_pow_n = (maxd*t)^n
//
//
#else
// Waterman method has difficulty with scale_pow_n overflowing for larger
// time values (~1-2)
//
// ...alternate method more robust to larger time values:
//
// for rate matrix Q:
// 1) maxd = -min(diag(Q))
// 2) set Q' = (Q + I.*maxd)
// 3) instead of exp(Qt) = sum_n(((Q.*t)^n)/(n!))
//    we calc    exp(Qt) = exp(-maxd*t) * sum_n(((Q'.*t)^n)/(n!))
//
// ..this is found from Q = Q' - I.*maxd
//                      exp(Q.*t) = exp(Q'.*t)*exp(-I.*maxd.*t)
//
// 4) we separate out time terms so that the time independent terms can
//   be cached for multi-time calls:
//               exp(Qt) = exp(-maxd*t) * sum_n((t^n)*(Q'^n)/(n!))
//
// note that in matrix_exp_scale below:
// emaxd = exp(-maxd*t)
// scale_pow_n = (t)^n
#endif
//
//
template <typename FloatType>
void
matrix_exp_series_prep(matrix_exp_series_prepdata<FloatType>& mepd,
                       const FloatType* in_mat){                 // dim*dim

  const unsigned N(mepd.dim);
  const unsigned N2(N*N);

  // huh... using symmetric matrix LAPACK solvers
  // seems to slow things down a bit:
  //mepd.is_symm = matrix_is_symm(in_mat,N);

  matrix_copy(in_mat,mepd.prep,N);

  // get maxd
  mepd.maxd = -mepd.prep[0];
  for(unsigned i(1);i<N;++i) {
    const FloatType v(-mepd.prep[i*(N+1)]);
    if(mepd.maxd<v) mepd.maxd=v;
  }

#ifdef ORIG_WATERMAN_METHOD
  // scale by 1/maxd
  //
  if(mepd.maxd != 0.) for(unsigned i(0);i<N2;++i) mepd.prep[i] /= mepd.maxd;
  // add I
  for(unsigned i(0);i<N;++i) mepd.prep[i*(N+1)] += 1.;
#else
  // add I.*maxd
  for(unsigned i(0);i<N;++i) mepd.prep[i*(N+1)] += mepd.maxd;
#endif

  mepd.info[0].is_calc = true;
  mepd.info[0].max = array_abs_max(mepd.prep,N2);
}




// get scaled matrix exponential from pre-calculated matrix taylor series.
//
template <typename FloatType>
void
matrix_exp_series_scale(FloatType* out_mat,    // dim*dim
                        matrix_exp_series_prepdata<FloatType>& mepd,
                        const FloatType in_scale,
                        workspace<char>& ws){

  const unsigned N = mepd.dim;
  const unsigned N2 = N*N;
#ifdef ORIG_WATERMAN_METHOD
  const FloatType maxd_scale(mepd.maxd*in_scale);
  const FloatType emaxd(std::exp(-maxd_scale));
#else
  const FloatType maxd_scale(in_scale);
  const FloatType emaxd(std::exp(-mepd.maxd*in_scale));
#endif

  if(emaxd<=0.|| is_float_invalid(emaxd)){
    std::ostringstream oss;
    oss << "matrix_exp_series_scale(): invalid_emaxd_value, in_scale, maxd_scale: "
        << emaxd << " " << in_scale << " " << maxd_scale;
    throw math_util_exception(oss.str().c_str());
  }

  matrix_identity(out_mat,N);

  ws.resize(mepd.size*sizeof(FloatType));
  FloatType* scale_pow(reinterpret_cast<FloatType* const>(ws.ptr()));

  // iter will be set dynamically according to precision requirements
  unsigned iter(0);
  FloatType last_iter_maxval(0);
  while(true){
    if(iter==0){
      scale_pow[iter] = maxd_scale;
    } else {
      scale_pow[iter] = scale_pow[iter-1]*maxd_scale;
      if( is_float_invalid(scale_pow[iter]) ){
        std::ostringstream oss;
        oss << "matrix_exp_series: invalid scale_pow value. iter/val/val[-1]/maxd_scale: "
            << iter << " " << scale_pow[iter] << " " << scale_pow[iter-1] << " " << maxd_scale;
        throw math_util_exception(oss.str().c_str());
      }
    }

    if( ! mepd.info[iter].is_calc ){
      scale_matrix_mult(static_cast<FloatType>(1./static_cast<FloatType>(iter+1)),
                        mepd.prep,
                        mepd.prep+((iter-1)*N2),
                        mepd.prep+(iter*N2),
                        N,
                        mepd.is_symm);

      mepd.info[iter].max = array_abs_max(mepd.prep+iter*N2,N2);

      if( is_float_invalid(mepd.info[iter].max) ){
        std::ostringstream oss;
        oss << "matrix_exp_series_scale(): invalid matrix max value. iter/max/max[-1]: "
            << iter << " " << mepd.info[iter].max << " " << mepd.info[iter-1].max;
        throw math_util_exception(oss.str().c_str());
      }
      mepd.info[iter].is_calc = true;
    }

    // maxval=mepd.max[iter]*scale_pow[iter]*emaxd is the maximum
    // difference in any component of the exp(Qt) approximation for
    // this iteration. Note that maxval for iteration i can be larger
    // than for iteration i-1, a condition which can occur before the
    // factorial in the denominator begins to dominate the iteration
    // (a transition referred to as 'the hump'). For this reason the
    // minimum threshold for maxval is not tested unless the value
    // of maxval for iteration i is less than that for i-1.
    //
#ifdef MEXP_DEBUG
    log_os << "iter,max,scale,max*scale,last max*scale,emaxd,in_scale,maxd: "
           << iter << " "
           << mepd.info[iter].max << " "
           << scale_pow[iter] << " "
           << mepd.info[iter].max*scale_pow[iter] << " "
           << last_iter_maxval << " "
           << emaxd << " "
           << in_scale << " "
           << mepd.maxd << "\n";
#endif

    const FloatType iter_maxval(mepd.info[iter].max*scale_pow[iter]);

    if(iter_maxval <= last_iter_maxval) {
      if(iter_maxval*emaxd<MATRIX_EXP_ABS_PREC) break;
    }

    last_iter_maxval = iter_maxval;

    iter++;
    if(mepd.size<=iter) {
      try {
        mepd.grow();
      } catch (math_util_exception& e){
        log_os << "MATH UTIL EXCEPTION: " << e.what() << "\n";
        log_os << "...caught in matrix_exp_series_scale(): iter_prec,goal_prec: "
               << iter_maxval*emaxd << " " << MATRIX_EXP_ABS_PREC << "\n";
        throw;
      }
      ws.resize_with_data(mepd.size*sizeof(FloatType));
      scale_pow = reinterpret_cast<FloatType* const>(ws.ptr());
    }
  }

  matrix_vector_mult_sum(mepd.prep,scale_pow,out_mat,N2,iter);

  for(unsigned i(0);i<N2;++i) out_mat[i] *= emaxd;
}

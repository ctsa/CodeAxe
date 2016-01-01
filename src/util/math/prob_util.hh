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

// $Id: prob_util.hh 871 2007-10-01 00:17:23Z ctsa $

/// \file

#include "math_util_exception.h"
#include "test_float.h"
#include "../general/metatags.h"

#include <cassert>
#include <cmath>

#include <iterator>
#include <limits>
#include <numeric>



template <typename ForwardIterator>
void
pdistro_norm(ForwardIterator pn,
             const ForwardIterator pn_end,
             const bool is_zero_safe){

  typedef typename std::iterator_traits<ForwardIterator>::value_type T;

  T sum(std::accumulate(pn,pn_end,T(0)));

  if(std::abs(sum) <= 0. && is_zero_safe){
    pdistro_unif(pn,pn_end);
    return;
  } else if(sum <= 0.){
    throw math_util_exception("pdistro_norm: invalid pdistro");
  }

  const T scale(T(1)/sum);
  for(;pn!=pn_end;++pn) *pn *= scale;
}



template <typename ForwardIterator1,
          typename ForwardIterator2>
void
pdistro_norm_free_param(ForwardIterator1 pn,
                        const ForwardIterator1 pn_end,
                        ForwardIterator2 fparam){

  {
    bool is_free(false);
    ForwardIterator2 fp(fparam);
    for(ForwardIterator1 p(pn);p!=pn_end;++p){
      if(*fp) { is_free=true; break; }
      ++fp;
    }
    if(! is_free) {
      pdistro_norm(pn,pn_end);
      return;
    }
  }

  typedef typename std::iterator_traits<ForwardIterator1>::value_type T;

  const T sum(std::accumulate(pn,pn_end,T(0)));
  if( sum < std::numeric_limits<T>::epsilon()){
    throw math_util_exception("pdistro_norm_free_param: invalid pdistro");
  }

  T freesum(0);
  ForwardIterator2 fp(fparam);
  for(ForwardIterator1 p(pn);p!=pn_end;++p){
    if(*fp) freesum = freesum + *p;
    ++fp;
  }
  if( freesum < std::numeric_limits<T>::epsilon()){
    throw math_util_exception("pdistro_norm_free_param: invalid free pdistro");
  }

  const T freegoal(T(1)-sum+freesum);
  if( freegoal < T(0) ){
    throw math_util_exception("pdistro_norm_free_param: invalid free pdistro");
  }

  const T scale(freegoal/freesum);
  for(;pn!=pn_end;++pn) {
    if(*fparam) *pn *= scale;
    ++fparam;
  }
}




void pdistro_check_fail_val(const double v,
                            const unsigned i)  NORETURN_TAG;

void pdistro_check_fail_sum(const double s) NORETURN_TAG ;


template <typename RandomAccessIterator>
void
pdistro_check(RandomAccessIterator pn,
              const unsigned N,
              double ptol,
              const unsigned stride){

  double sum(0.);
  for(unsigned i(0);i<N;++i){
    const unsigned idx = i*stride;
    const double val(*(pn+idx));
    if(val<0.|| is_float_invalid(val)){
      pdistro_check_fail_val(val,idx);
    }
    sum += val;
  }

  if(std::fabs(sum-1.)>ptol){
    pdistro_check_fail_sum(sum);
  }
}



// T must be some floating-point type
template <typename ForwardIterator>
void
pdistro_unif(ForwardIterator pn,
             const ForwardIterator pn_end){
  typedef typename std::iterator_traits<ForwardIterator>::value_type T;

  const T val(T(1)/static_cast<T>(std::distance(pn,pn_end)));
  std::fill(pn,pn_end,val);
}

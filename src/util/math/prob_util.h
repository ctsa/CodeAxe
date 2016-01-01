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

// $Id: prob_util.h 1198 2008-04-21 23:09:14Z ctsa $

/// \file
///
/// \brief simple discrete probability distro utilities
///


#ifndef __PROB_UTIL_H
#define __PROB_UTIL_H

#include <cmath>

#include <algorithm>
#include <iterator>


/// \brief normalize distro by adjusting all parameters
///
/// if is_zero_safe, then a uniform distro is returned
/// when the input distro sums to a nonpositive value
///
template <typename ForwardIterator>
void
pdistro_norm(ForwardIterator pn,
             const ForwardIterator pn_end,
             const bool is_zero_safe=false);

template <typename ForwardIterator>
void
pdistro_norm_safe(ForwardIterator pn,
                  const ForwardIterator pn_end){
  pdistro_norm(pn,pn_end,true);
}

/// normalize distribution by only adjusting free parameters. if
/// all parameters are locked, adjust all. returns scaling factor
///
template <typename ForwardIterator1,
          typename ForwardIterator2>
void
pdistro_norm_free_param(ForwardIterator1 pn,
                        const ForwardIterator1 pn_end,
                        ForwardIterator2 fparam);

/// check that a probability distro sums to unity
/// at the specified tolerance
///
template <typename RandomAccessIterator>
void
pdistro_check(RandomAccessIterator pn,
              const unsigned N,
              double ptol,
              const unsigned stride=1);

template <typename ForwardIterator>
void
pdistro_unif(ForwardIterator pn,
             const ForwardIterator pn_end);

/// safe for inplace case (pdf == cdf)
template <typename FloatType>
void
pdistro_to_cdf(const FloatType pdf[],
               FloatType cdf[],
               const unsigned N){
  cdf[0] = pdf[0];
  for(unsigned i(1);i<N;++i) cdf[i] = cdf[i-1]+pdf[i];
}

#if 0
template <typename ForwardIterator>
void
pdistro_to_cdf_inplace(ForwardIterator p,
                       const ForwardIterator p_end){
  if(p!=p_end){
    ForwardIterator plast(p);
    while((++p)!=p_end){
      *p += *plast;
      plast = p;
    }
  }
}
#endif

template <typename RandomAccessIterator>
void
pdistro_to_cdf_inplace(RandomAccessIterator p,
                       const unsigned N){
  for(unsigned i(1);i<N;++i) *(p+i) += *(p+i-1);
}


/// safe for inplace case (cdf == pdf)
template <typename FloatType>
void
pdistro_from_cdf(const FloatType cdf[],
                 FloatType pdf[],
                 const unsigned N){
  for(unsigned i(N-1);i>0;--i) pdf[i] = cdf[i]-cdf[i-1];
  pdf[0] = cdf[0];
}


template <typename RandomAccessIterator>
void
pdistro_from_cdf_inplace(RandomAccessIterator p,
                         const unsigned N){
  for(unsigned i(N-1);i>0;--i) *(p+i) -= *(p+i-1);
}


/// relative entropy of two distros:
template <typename FloatType>
FloatType
rel_ent(const FloatType* pdf1,
        const FloatType* pdf2,
        const unsigned N){

  FloatType re(0.);
  for(unsigned i(0);i<N;++i){
    if(fabs(pdf1[i]) <= 0.) continue; // 0*-inf=0
    re += pdf1[i]*std::log(pdf1[i]/pdf2[i]);
  }
  return re;
}


#include "prob_util.hh"

#endif

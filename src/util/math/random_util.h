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

// $Id: random_util.h 1148 2008-03-11 00:49:49Z ctsa $

/// \file

#ifndef __RANDOM_UTIL_H
#define __RANDOM_UTIL_H

#define USE_BOOST_RANDOM

#ifdef USE_BOOST_RANDOM
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#else
#include <cstdlib>
#endif

#include <cassert>
#include <cmath>

#include <algorithm>
#include <vector>

#ifdef USE_BOOST_RANDOM
typedef boost::mt19937 base_generator_type;
typedef boost::variate_generator<base_generator_type&, boost::uniform_real<> > br_uni_func_type;
extern br_uni_func_type br_uni_func;
#endif

long int
random_init(long int seed = -1);

inline
double
random_uniform(){
#ifdef USE_BOOST_RANDOM
  return br_uni_func();
#else
  return drand48();
#endif
}

inline
double
random_exponential(double l){
  double r;
  while((r=random_uniform())<=0.){};
  return -std::log(r)/l;
}

/// \brief pull a random variate from a discrete cdf
///
template <typename FloatType>
unsigned
random_cdf_variate(const FloatType cdf[],
                   const unsigned N,
                   const FloatType scale=1.,
                   const FloatType offset=0.){

  const FloatType r(static_cast<FloatType>(offset+(scale*random_uniform())));
  const FloatType* lbp(std::lower_bound(cdf,cdf+N,r));
  return std::min(static_cast<unsigned>(lbp-cdf),N-1);
}


/// \brief pull a random variate from an unnormalized discrete cdf
///
template <typename FloatType>
unsigned
random_scaled_cdf_variate(const FloatType scdf[],
                          const unsigned N){

  return random_cdf_variate(scdf,N,scdf[N-1]);
}


/// \brief pull a random variate from a discrete cdf in [lower,upper)
///
template <typename FloatType>
unsigned
random_subset_cdf_variate(const FloatType cdf[],
                          const unsigned lower,
                          const unsigned upper){
  assert(lower<upper);

  FloatType offset(0);
  if(lower>1) offset=cdf[lower-1];

  const FloatType scale(cdf[upper-1]-offset);

  return lower+random_cdf_variate(cdf+lower,upper-lower,scale,offset);
}


/// \brief a discrete distribution which allows efficient updates to
/// single state probs
///
/// efficiently produce random variates from a very large unnormalized
/// probability distribution subject to frequent single state value
/// changes by breaking the equiv cdfs into a tree structure. This allows
/// [?? check this later->] O(N log(N)) updates
///
template <typename FloatType>
struct random_scaled_pdistro_variate_cached {

  random_scaled_pdistro_variate_cached(const FloatType* scaled_pdistro,
                                       const unsigned N,
                                       const unsigned max_cdf_size=1024);

  ~random_scaled_pdistro_variate_cached();

  unsigned get_random_draw() const {
    unsigned combined_draw(0);

    for(unsigned i(0);i<_n_cdf_draws;++i){
      const unsigned cdf_start(combined_draw ? combined_draw/_draw_denom[i-1] : 0);
      FloatType* reduced_cdf(_cache[i]+cdf_start*_max_cdf_size);
      const unsigned reduced_size(std::min(div_ceil(_size-combined_draw,_draw_denom[i]),_max_cdf_size));

      const unsigned draw(random_scaled_cdf_variate(reduced_cdf,reduced_size));

      combined_draw += draw*_draw_denom[i];
    }
    assert(combined_draw<_size);
    return combined_draw;
  }

  void update(const unsigned pos,
              const FloatType delta) {
    unsigned cell_start(0);

    for(unsigned i(0);i<_n_cdf_draws;++i){
      const unsigned cdf_start(cell_start ? cell_start/_draw_denom[i-1] : 0);
      FloatType* reduced_cdf(_cache[i]+cdf_start*_max_cdf_size);
      const unsigned reduced_size(std::min(div_ceil(_size-cell_start,_draw_denom[i]),_max_cdf_size));

      const unsigned cell((pos-cell_start)/_draw_denom[i]);
      for(unsigned j(cell);j<reduced_size;++j) { reduced_cdf[j] += delta; }
      cell_start += cell*_draw_denom[i];
    }
  }

private:
  // returns equivalent of int(ceil(float(x)/float(y))) ...or same as
  // '/>' operator in scheme
  //
  static unsigned div_ceil(const unsigned x,
                           const unsigned y){
    if(x) return 1+((x-1)/y);
    else  return 0;
  }

  unsigned _size;
  unsigned _max_cdf_size;
  unsigned _n_cdf_draws;
  std::vector<unsigned> _draw_denom;
  FloatType** _cache;
};


#include "random_util.hh"

#endif

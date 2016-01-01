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

// $Id: random_util.hh 743 2007-08-14 15:47:12Z ctsa $

/// \file

template <typename FloatType>
random_scaled_pdistro_variate_cached<FloatType>::
random_scaled_pdistro_variate_cached(const FloatType* scaled_pdistro,
                                     const unsigned N,
                                     const unsigned max_cdf_size)
  : _size(N), _max_cdf_size(std::min(_size,max_cdf_size)),
    _n_cdf_draws(0), _draw_denom(), _cache(0) {

  assert(N>0);
  assert(max_cdf_size>1);

  for(unsigned i(_size);i>1;i=div_ceil(i,_max_cdf_size)){ _n_cdf_draws++; }

  _cache = new FloatType*[_n_cdf_draws];
  _draw_denom.resize(_n_cdf_draws);

  for(unsigned i(0);i<_n_cdf_draws;++i){
    unsigned& draw_denom(_draw_denom[i]);
    draw_denom=1;
    {
      const unsigned draw_denom_pow((_n_cdf_draws-1-i));
      for(unsigned j(0);j<draw_denom_pow;++j) { draw_denom *= _max_cdf_size; }
    }

    const unsigned n_cdfs(i ? div_ceil(_size,_draw_denom[i-1]) : 1);

    _cache[i] = new FloatType[n_cdfs*_max_cdf_size];

    for(unsigned c(0);c<n_cdfs;++c){
      FloatType* reduced_cdf(_cache[i]+c*_max_cdf_size);

      const unsigned fd(c ? c*_draw_denom[i-1] : 0);

      const unsigned reduced_size( (c+1)==n_cdfs
                                   ? div_ceil(_size-fd,draw_denom)
                                   : _max_cdf_size);

      if(draw_denom>1){
        for(unsigned j(0);j<reduced_size;++j){
          if(j==0) reduced_cdf[j]=0.;
          else     reduced_cdf[j]=reduced_cdf[j-1];

          const unsigned joffset(fd+j*draw_denom);
          const unsigned kend(std::min(_size-joffset,draw_denom));
          for(unsigned k(0);k<kend;++k){
            reduced_cdf[j]+=scaled_pdistro[joffset+k];
          }
        }
      } else {
        pdistro_to_cdf(scaled_pdistro+fd,reduced_cdf,reduced_size);
      }
    }
  }
}



template <typename FloatType>
random_scaled_pdistro_variate_cached<FloatType>::
~random_scaled_pdistro_variate_cached(){
  if(_cache){
    for(unsigned i(0);i<_n_cdf_draws;++i){
      if(_cache[i]) delete [] _cache[i];
    }
    delete [] _cache;
  }
}

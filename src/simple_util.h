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

// $Id: simple_util.h 1148 2008-03-11 00:49:49Z ctsa $

/// \file

#ifndef __SIMPLE_UTIL_H
#define __SIMPLE_UTIL_H

#include <algorithm>

inline
unsigned
off_diag_size(const unsigned a){
  return (a*(a-1));
}

#if 0
inline
unsigned
encode_nmer(const unsigned a[],
            const unsigned len,
            const unsigned N){
  unsigned x = a[0];
  for(unsigned i(1);i<len;++i) x = x*N+a[i];
  return x;
}


inline
void
decode_nmer(const unsigned site_id,
            const unsigned len,
            const unsigned N,
            unsigned a[]){
  unsigned x = site_id;
  for(unsigned i(len-1);i;--i){
    a[i] = x%N;
    x /= N;
  }
  a[0] = x;
}
#endif


inline
unsigned
periodic_decrement(const unsigned val,
                   const unsigned period){
  //assert(period>0);
  return (val+(period-1))%period;
}


inline
unsigned
periodic_increment(const unsigned val,
                   const unsigned period){
  //assert(period>0);
  return (val+1)%period;
}



template <typename T>
struct simple_array_base {
protected:
  simple_array_base(const unsigned d1=0) : val(0), _dim1(d1) {}

  simple_array_base(const simple_array_base<T>& s)
    : val(0), _dim1(s._dim1) {
    alloc();
    for(unsigned i(0);i<_dim1;++i) val[i] = s.val[i];
  }

  ~simple_array_base(){ dealloc(); }

public:
  unsigned dim1() const { return _dim1; }

  T& operator[](const unsigned i) { return val[i]; }
  const T& operator[](const unsigned i) const { return val[i]; }

  void set_val(const T v){ std::fill(val,val+size(),v); }

  T* ptr() { return val; }
  const T * ptr() const { return val; }

  unsigned size() const { return dim1(); }

  T* begin() { return val; }
  const T * begin() const { return val; }
  T* end() { return val+size(); }
  const T * end() const { return val+size(); }

  void clear() { dealloc(); _dim1=0; }

protected:
  void alloc(){ if(_dim1) val = new T[_dim1]; }

  void dealloc(){ if(_dim1) {delete [] val; val=0;} }

  T* val;
  unsigned _dim1;
};


template <typename T>
struct simple_array : public simple_array_base<T> {
  typedef simple_array_base<T> base_t;

  simple_array(const unsigned d1) : base_t(d1) { base_t::alloc(); }

  simple_array(const unsigned d1,
               const T init_val) : base_t(d1) {
    base_t::alloc();
    set_val(init_val);
  }

private:
  simple_array& operator=(const simple_array<T>&);
};



template <typename T>
struct simple_init_array : public simple_array_base<T> {
  typedef simple_array_base<T> base_t;

  simple_init_array& operator=(const simple_init_array<T>& rhs) {
    if(&rhs != this) deep_copy(rhs);
    return *this;
  }

  void init(const unsigned d1){
    if(base_t::val && base_t::_dim1 == d1) return;
    base_t::dealloc();
    base_t::_dim1 = d1;
    base_t::alloc();
  }

  void init(const unsigned d1,
            const T init_val){
    init(d1);
    base_t::set_val(init_val);
  }

private:
  void deep_copy(const simple_init_array<T>& s){
    init(s._dim1);
    for(unsigned i(0);i<base_t::_dim1;++i) base_t::val[i] = s.val[i];
  }
};


template <typename T>
struct simple_matrix_base {
protected:
  simple_matrix_base(const unsigned d1=0,
                     const unsigned d2=0)
    : val(0), _dim1(d1), _dim2(d2) {}

  simple_matrix_base(const simple_matrix_base<T>& s)
    : val(0), _dim1(s._dim1), _dim2(s._dim2) {
    alloc();
    for(unsigned i(0);i<_dim1*_dim2;++i) val[0][i] = s.val[0][i];
  }

  ~simple_matrix_base(){ dealloc(); }

public:
  unsigned dim1() const { return _dim1; }
  unsigned dim2() const { return _dim2; }

  T*& operator[](const unsigned i) { return val[i]; }
  const T * operator[](const unsigned i) const { return val[i]; }

  void set_val(const T v){ if(_dim1) std::fill(*val,*val+size(),v); }

  T** ptr() { return val; }
  const T * const * ptr() const { return val; }

  void clear() { dealloc(); _dim1=0; _dim2=0; }

private:
  unsigned size() const { return dim1()*dim2(); }

protected:
  void alloc(){
    if(_dim1) {
      val = new T*[_dim1];
      if(_dim2) {
        val[0] = new T[_dim1*_dim2];
      }
      for(unsigned i(1);i<_dim1;++i) val[i] = val[0]+i*_dim2;
    }
  }

  void dealloc() {
    if(_dim1){
      if(_dim2){
        delete [] val[0];
      }
      delete [] val;
      val=0;
    }
  }

  T** val;
  unsigned _dim1;
  unsigned _dim2;
};


template <typename T>
struct simple_matrix : public simple_matrix_base<T> {
  typedef simple_matrix_base<T> base_t;

  simple_matrix(const unsigned d1,
                const unsigned d2)
    : base_t(d1,d2) { base_t::alloc(); }

  simple_matrix(const unsigned d1,
                const unsigned d2,
                const T init_val)
    : base_t(d1,d2) { base_t::alloc(); set_val(init_val); }

private:
  simple_matrix& operator=(const simple_matrix<T>&);
};



template <typename T>
struct simple_init_matrix : simple_matrix_base<T> {
  typedef simple_matrix_base<T> base_t;

  simple_init_matrix& operator=(const simple_init_matrix<T>& rhs) {
    if(&rhs != this) deep_copy(rhs);
    return *this;
  }

  void init(const unsigned d1,
            const unsigned d2) {
    if(base_t::val && base_t::_dim1 == d1 && base_t::_dim2 == d2) return;
    base_t::dealloc();
    base_t::_dim1 = d1;
    base_t::_dim2 = d2;
    base_t::alloc();
  }

  void init(const unsigned d1,
            const unsigned d2,
            const T init_val) {
    init(d1,d2);
    set_val(init_val);
  }

private:
  void deep_copy(const simple_init_matrix<T>& s){
    init(s._dim1,s._dim2);
    for(unsigned i(0);i<base_t::_dim1*base_t::_dim2;++i) base_t::val[0][i] = s.val[0][i];
  }
};



template <typename T>
struct simple_matrix3d_base {
protected:
  simple_matrix3d_base(const unsigned d1 = 0,
                       const unsigned d2 = 0,
                       const unsigned d3 = 0)
  : val(0), _dim1(d1), _dim2(d2), _dim3(d3) {}

  simple_matrix3d_base(const simple_matrix3d_base<T>& s)
    : val(0), _dim1(s._dim1), _dim2(s._dim2), _dim3(s._dim3) {
    alloc();
    for(unsigned i(0);i<_dim1*_dim2*_dim3;++i) val[0][0][i] = s.val[0][0][i];
  }

  ~simple_matrix3d_base(){ dealloc(); }

private:
  simple_matrix3d_base& operator=(const simple_matrix3d_base<T>&);

public:
  unsigned dim1() const { return _dim1; }
  unsigned dim2() const { return _dim2; }
  unsigned dim3() const { return _dim3; }

  T**& operator[](const unsigned i) { return val[i]; }
  const T* const * operator[](const unsigned i) const { return val[i]; }

  void set_val(const T v){ if(_dim1 && _dim2) std::fill(**val,**val+size(),v); }

  T*** ptr() { return val; }
  const T * const * const * ptr() const { return val; }

  void clear() { dealloc(); _dim1=0; _dim2=0; _dim3=0; }

private:
  unsigned size() const { return dim1()*dim2()*dim3(); }

protected:
  void alloc(){
    if(_dim1){
      val = new T**[_dim1];
      if(_dim2){
        val[0] = new T*[_dim1*_dim2];
        if(_dim3){
          val[0][0] = new T[_dim1*_dim2*_dim3];
        }
      }
    }
    for(unsigned i(0);i<_dim1;++i){
      val[i] = val[0]+i*_dim2;
      for(unsigned j(0);j<_dim2;++j) {
        val[i][j] = val[0][0] + (i*_dim2+j)*_dim3;
      }
    }
  }

  void dealloc(){
    if(_dim1){
      if(_dim2) {
        if(_dim3) {
          delete [] val[0][0];
        }
        delete [] val[0];
      }
      delete [] val;
      val=0;
    }
  }

  T*** val;
  unsigned _dim1;
  unsigned _dim2;
  unsigned _dim3;
};



template <typename T>
struct simple_matrix3d : public simple_matrix3d_base<T> {

  typedef simple_matrix3d_base<T> base_t;

  simple_matrix3d(const unsigned d1,
                  const unsigned d2,
                  const unsigned d3)
    : base_t(d1,d2,d3) { base_t::alloc(); }

  simple_matrix3d(const unsigned d1,
                  const unsigned d2,
                  const unsigned d3,
                  const T init_val)
    : base_t(d1,d2,d3) { base_t::alloc(); set_val(init_val); }
};



template <typename T>
struct simple_init_matrix3d : public simple_matrix3d_base<T> {

  typedef simple_matrix3d_base<T> base_t;

  void init(const unsigned d1,
            const unsigned d2,
            const unsigned d3){
    if(base_t::val && base_t::_dim1 == d1 && base_t::_dim2 == d2 && base_t::_dim3 == d3) return;
    base_t::dealloc();
    base_t::_dim1 = d1;
    base_t::_dim2 = d2;
    base_t::_dim3 = d3;
    base_t::alloc();
  }

  void init(const unsigned d1,
            const unsigned d2,
            const unsigned d3,
            const T init_val){
    init(d1,d2,d3);
    set_val(init_val);
  }
};

#endif

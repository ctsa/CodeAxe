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

// $Id: indy_random_var.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __INDY_RANDOM_VAR_H
#define __INDY_RANDOM_VAR_H


/// \brief a class to automatically track the variance and mean of independent
/// random variables through simple mathematical ops
///
template <typename T>
struct irv_t {
  typedef irv_t<T> self_t;

  irv_t(T init_m=0,
        T init_v=0) : m(init_m), v(init_v) {}


  const self_t operator+(const self_t& rhs) const {
    return self_t(m+rhs.m,v+rhs.v);
  }

  const self_t operator-(const self_t& rhs) const {
    return self_t(m-rhs.m,v+rhs.v);
  }

  const self_t operator*(const self_t& rhs) const {
    // based on first order taylor approx..
    return self_t(m*rhs.m,(rhs.m*rhs.m)*v+(m*m)*rhs.v);
  }

  const self_t operator/(const self_t& rhs) const {
    // based on first order taylor approx..
    return self_t(m/rhs.m,v/(rhs.m*rhs.m)+(m*m)*rhs.v/(rhs.m*rhs.m*rhs.m*rhs.m));
  }

  const self_t operator+=(const self_t& rhs){
    *this = *this+rhs;
    return *this;
  }

  const self_t operator-=(const self_t& rhs){
    *this = *this-rhs;
    return *this;
  }

  const self_t operator*=(const self_t& rhs){
    *this = *this*rhs;
    return *this;
  }

  const self_t operator/=(const self_t& rhs){
    *this = *this/rhs;
    return *this;
  }

  operator T() const { return m; }

  T m;
  T v;
};


// scalar interaction
template <typename T>
const irv_t<T> operator+(const irv_t<T>& lhs,
                         const T& rhs){
  return irv_t<T>(lhs.m*rhs,lhs.v);
}

template <typename T>
const irv_t<T> operator+(const T& lhs,
                         const irv_t<T>& rhs){
  return rhs+lhs;
}

template <typename T>
const irv_t<T> operator*(const irv_t<T>& lhs,
                         const T& rhs){
  return irv_t<T>(lhs.m*rhs,lhs.v*rhs*rhs);
}

template <typename T>
const irv_t<T> operator*(const T& lhs,
                         const irv_t<T>& rhs){
  return rhs*lhs;
}

template <typename T>
const irv_t<T> operator/(const irv_t<T>& lhs,
                         const T& rhs){
  return lhs*(1/rhs);
}

#endif

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// SubsTK : phylogenetic analysis and simulation library
//
//   http://www.phrap.org
//
//
// Copyright (C) 2007 Christopher T Saunders (ctsa@u.washington.edu)
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

// $Id: workspace.h 1055 2007-12-08 04:27:31Z ctsa $

/// \file

#ifndef __WORKSPACE_H
#define __WORKSPACE_H


#include "uncopyable.h"


template <typename T>
struct workspace : private uncopyable {

  explicit
  workspace(const unsigned n = 0) : _ptr(0), _size(0) {
    if(n) resize(n);
  }

  ~workspace() { if(_ptr) delete [] _ptr; }

  void
  resize(const unsigned s){
    if(_size<s){
      _size=s;
      if(_ptr) delete [] _ptr;
      _ptr = new T[_size];
    }
  }

  void
  resize_with_data(const unsigned s){
    if(_size<s){
      T* new_ptr(new T[s]);
      std::copy(_ptr,_ptr+_size,new_ptr);
      _size=s;
      if(_ptr) delete [] _ptr;
      _ptr = new_ptr;
    }
  }

  T* ptr() { return _ptr; }

private:
  T* _ptr;
  unsigned _size;
};

#endif

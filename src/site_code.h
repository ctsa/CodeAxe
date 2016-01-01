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

// $Id: site_code.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __SITE_CODE_H
#define __SITE_CODE_H

#include "util/general/die.h"

#include <cassert>
#include <cstdlib>

#include <iosfwd>


struct site_code {
  typedef site_code self_t;
  typedef unsigned short index_type;

  site_code(const index_type ntaxa=0)
    : _ntaxa(ntaxa), _site(0) { init(); }

  site_code(const self_t& s)
    : _ntaxa(s._ntaxa), _site(0) {
    init();
    for(index_type i(0);i<_ntaxa;++i){ _site[i] = s._site[i]; }
  }

  self_t& operator=(const self_t& rhs){
    if(&rhs==this) return *this;
    set_ntaxa(rhs._ntaxa);
    for(index_type i(0);i<_ntaxa;++i){ _site[i] = rhs._site[i]; }
    return *this;
  }

  ~site_code() { clear(); }

  void
  set_ntaxa(const index_type ntaxa=0){
    if( ntaxa != _ntaxa ){
      clear();
      _ntaxa = ntaxa;
      init();
    }
  }

  index_type get_ntaxa() const { return _ntaxa; }

  void
  set_taxid(const index_type t, const index_type site_val){
    assert(t<_ntaxa);
    _site[t]=site_val;
  }

  index_type
  get_taxid(const index_type t) const {
    assert(t<_ntaxa);
    return _site[t];
  }

  void
  check_init(const unsigned n_seqs,
             const index_type s) {
    if(get_ntaxa()==0){
      set_ntaxa(n_seqs);
      for(unsigned i(0);i<n_seqs;++i) set_taxid(i,s);
    }
  }


  bool
  operator==(const site_code& right) const {
    if(_ntaxa != right._ntaxa) return false;
    for(index_type i(0);i<_ntaxa;++i){
      if( _site[i] != right._site[i]) return false;
    }
    return true;
  }

  bool
  operator<(const site_code& right) const {
    if(_ntaxa != right._ntaxa) { die("attempting to compare site_codes with different numbers of sequences"); };
    for(index_type i(0);i<_ntaxa;++i){
      if( _site[i] != right._site[i]) return  _site[i] < right._site[i];
    }
    return false;
  }

private:
  void clear() {  if(_site) delete [] _site; _site=0; }
  void init() { if(_ntaxa) _site = new index_type[_ntaxa]; }

  index_type _ntaxa;
  index_type* _site;
};

std::ostream& operator<<(std::ostream& os,const site_code& sc);







struct cat_site_code {

  typedef unsigned char cat_type;

  bool
  operator<(const cat_site_code& right) const {
    if(site == right.site){
      return cat_no < right.cat_no;
    } else {
      return site < right.site;
    }
  }

  site_code site;
  cat_type cat_no;
};







#if 0
struct multi_site {
  multi_site() {}

  multi_site(const site_code& _prev,
             const site_code& _current)
    : prev(_prev), current(_current) {}

  bool
  operator<(const multi_site& right) const {
    if(current == right.current) {
      return prev < right.prev;
    } else {
      return current < right.current;
    }
  }

  site_code prev;
  site_code current;
};
#endif




#endif

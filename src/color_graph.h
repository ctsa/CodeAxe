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

// $Id: color_graph.h 1037 2007-12-03 23:58:34Z ctsa $

/// \file

#ifndef __COLOR_GRAPH_H
#define __COLOR_GRAPH_H

#include "util/general/uncopyable.h"

#include <cassert>

#include <map>
#include <vector>



struct color_vertex : private uncopyable {

  color_vertex() : color(0), is_colored(false) {}

  void
  add_edge(color_vertex* to) { _outedge.push_back(to); }

  void
  color_connected(const unsigned c){
    if(! is_colored){
      is_colored=true;
      color=c;
      const unsigned es(_outedge.size());
      for(unsigned i(0);i<es;++i) _outedge[i]->color_connected(c);
    }
  }

  unsigned color;
  bool is_colored;
private:
  std::vector<color_vertex*> _outedge;
};



template <typename I>
class color_graph : private uncopyable {
private:
  typedef std::map<I,color_vertex*> vlup_t;
public:
  color_graph() {}

  ~color_graph() {
    typename vlup_t::iterator i=_vlup.begin(),i_end=_vlup.end();
    for(;i!=i_end;++i) delete i->second;
  }

  // nothing prevents repeated edges:
  void
  add_edge(const I& i0,const I& i1){
    assert(i0 != i1);

    color_vertex* v0(get_vertex(i0));
    color_vertex* v1(get_vertex(i1));

    v0->add_edge(v1);
    v1->add_edge(v0);
  }

  bool
  vertex_is_colored(const I& i){ return get_vertex(i)->is_colored; }

  unsigned
  vertex_color(const I& i){ return get_vertex(i)->color; }

  void
  color_connected(const I& i,const unsigned c){
    get_vertex(i)->color_connected(c);
  }

private:
  color_vertex*
  get_vertex(const I& i){
    const typename vlup_t::const_iterator a(_vlup.find(i));
    if( a == _vlup.end() ){
      _vlup[i] = new color_vertex();
      return _vlup[i];
    }
    return a->second;
  }

  vlup_t _vlup;
};


#endif

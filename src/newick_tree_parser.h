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

// $Id: newick_tree_parser.h 1184 2008-03-27 21:15:54Z ctsa $

/// \file

#ifndef __NEWICK_TREE_PARSER_H
#define __NEWICK_TREE_PARSER_H

#include "subs_ml_types.h"
#include "substk_exception.h"

#include "util/general/uncopyable.h"

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>


struct newick_tree_node;


/// \brief parses newick tree string into a simple tree structure
///
/// note that spaces (but not other w/s) are allowed within the branch
/// label
///
/// \todo revamp error handling to return error codes rather then
/// throw exceptions, b/c invalid tree input is not an unexpected event
///
struct newick_tree_parser : private uncopyable {

  newick_tree_parser(const char* tree_string);

  const newick_tree_node* root() const { return _root.get(); }

private:
  std::auto_ptr<newick_tree_node> _root;
};


std::ostream& operator<<(std::ostream& os,
                         const newick_tree_parser& b);



struct newick_tree_node {

  newick_tree_node() : value(0.), is_value_set(false) {}

  ~newick_tree_node() {
    const unsigned csize(_child.size());
    for(unsigned i(0);i<csize;++i) delete _child[i];
  }

private:
  newick_tree_node(const newick_tree_node&);
  newick_tree_node& operator=(const newick_tree_node&);
public:

  unsigned child_size() const { return _child.size(); }

  const newick_tree_node* child(const unsigned i) const { return _child[i]; }

  void add_child(newick_tree_node* b){
    static const char* const msg("newick_tree_node.add_child(): invalid argument");

    if(b==0 || b==this) throw substk_exception(msg);

    const unsigned cs(child_size());
    for(unsigned i(0);i<cs;++i) if(b==child(i)) throw substk_exception(msg);

    _child.push_back(b);
  }

  std::string label;
  smlfloat value;
  bool is_value_set;

private:
  std::vector<newick_tree_node*> _child;
};

#endif

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

// $Id: bi_tree.h 999 2007-11-02 18:48:00Z ctsa $

/// \file

#ifndef __BI_TREE_H
#define __BI_TREE_H

#include "substk_exception.h"

#include <cassert>

#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <utility>


struct bi_tree_node;



/// \breif indexed bifurcating tree
///
/// \todo revamp error handling in the string ctor to return error
/// codes rather then throw exceptions, b/c invalid tree input is not
/// unexpected
///
struct bi_tree {

  bi_tree(const char* newick_tree_string);

  bi_tree(const bi_tree& s);

private:
  bi_tree& operator=(const bi_tree&);
public:

  unsigned leaf_size() const { return _leaf_node_index.size(); }
  unsigned branch_size() const { return ( node_size()==0 ? 0 : node_size()-1); }
  unsigned node_size() const { return _node_index.size(); }

  const bi_tree_node* root() const { return _root.get(); }

  const bi_tree_node* node(const std::string& s) const {
    std::map<std::string,const bi_tree_node*>::const_iterator i(_node_label_map.find(s));
    if(i==_node_label_map.end()){
      throw substk_exception("bi_tree.node(std::string): invalid node name");
    }
    return i->second;
  }

  bool is_node_label(const std::string& s) const {
    return ! ( _node_label_map.find(s) == _node_label_map.end());
  }

  const bi_tree_node* node(const unsigned i) const { return _node_index[i]; }
  const bi_tree_node* branch_node(const unsigned i) const { return node(i+1); }
  const bi_tree_node* leaf_node(const unsigned i) const { return _leaf_node_index[i]; }

  bi_tree_node* branch_node(const unsigned i) { return _node_index[i+1]; }
  void clear_branch_lengths(const bool is_warn = true);

private:
  bi_tree_node* root() { return _root.get(); }

  void ctor_setup_tree_index();


  std::auto_ptr<bi_tree_node> _root;
  std::vector<bi_tree_node*> _node_index;
  std::vector<const bi_tree_node*> _leaf_node_index;
  std::map<std::string,const bi_tree_node*> _node_label_map;
};



std::ostream& operator<<(std::ostream& os,
                         const bi_tree& b);



struct newick_tree_node;

/// bifurcating tree node, with redundant parent edges included
/// to make things easier
///
struct bi_tree_node {

private:
  friend struct bi_tree;

  /// both ctors copy the entire tree structure below them:
  ///
  bi_tree_node(const newick_tree_node* b);
  bi_tree_node(const bi_tree_node* b);

  bi_tree_node(const bi_tree_node&);
  bi_tree_node& operator=(const bi_tree_node&);

public:

  void set_length(const double l) {
    _data.length=l;
    _data.is_length_set=true;
  }

  const std::string& label() const { return _data.label; }
  double length() const { return _data.length; }
  bool is_length_set() const { return _data.is_length_set; }
  short int node_id() const { return _data.node_id; }
  short int branch_id() const {
    if(node_id()<=0) die("invalid call ti bi_tree_node.branch_id()");
    return node_id()-1;
  }
  short int leaf_id() const { return _data.leaf_id; }

  const bi_tree_node* parent() const { return _parent; }
  const bi_tree_node* child1() const { return _child1.get(); }
  const bi_tree_node* child2() const { return _child2.get(); }

  const bi_tree_node* sister() const {
    if     (_parent==0)              return 0;
    else if(_parent->child2()==this) return _parent->child1();
    else                             return _parent->child2();
  }

  bool is_root() const { return parent()==0; }
  bool is_leaf() const { return child1()==0; }

private:
  bi_tree_node* child1() { return _child1.get(); }
  bi_tree_node* child2() { return _child2.get(); }

  void clear_length(const bool is_warn);

  void make_index(std::vector<bi_tree_node*>& node_index,
                  std::vector<const bi_tree_node*>& leaf_node_index);

  void set_children(bi_tree_node* b1,
                    bi_tree_node* b2){
    if(b1==0 || b2==0 || b1 == this || b2 == this || b1 == b2) {
      die("bi_tree_node.set_children(): invalid arguments");
    }

    _child1.reset(b1); child1()->_parent=this;
    _child2.reset(b2); child2()->_parent=this;
  }

  struct auto_copy {
    auto_copy()
      : length(0.),is_length_set(false),
        node_id(-1), leaf_id(-1) {}

    std::string label;
    double length;
    bool is_length_set;
    short int node_id;
    short int leaf_id;
  };

  auto_copy _data;
  bi_tree_node* _parent;
  std::auto_ptr<bi_tree_node> _child1;
  std::auto_ptr<bi_tree_node> _child2;
};

#endif

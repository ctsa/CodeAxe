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

// $Id: bi_tree.cc 984 2007-10-31 22:00:31Z ctsa $

/// \file

#include "bi_tree.h"
#include "newick_tree_parser.h"

#include "util/general/die.h"

#include <ostream>
#include <sstream>



bi_tree_node::
bi_tree_node(const newick_tree_node* b)
  : _data(), _parent(0) {
  const unsigned nc(b->child_size());
  if(! ( nc == 0 || nc == 2 ) ) {
    throw substk_exception("bi_tee_node(): bi_tree_type requires a bifurcating tree");
  }

  if(nc){
    set_children(new bi_tree_node(b->child(0)),new bi_tree_node(b->child(1)));
  }

  if(b->label.find_first_of(" ") != std::string::npos){
    throw substk_exception("bi_tree_node(): bi_tree_type does not accept w/s in tree node labels");
  }
  _data.label=b->label;
  _data.is_length_set=b->is_value_set;
  _data.length=b->value;
}



bi_tree_node::
bi_tree_node(const bi_tree_node* b)
  : _data(b->_data), _parent(0) {

  if(b==0) die("bi_tree_node(): invalid ctor argument");
  if(b->is_leaf()) return;
  set_children(new bi_tree_node(b->child1()),new bi_tree_node(b->child2()));
}



void
bi_tree_node::
clear_length(const bool is_warn){

  if(is_length_set()) {
    if(is_warn) warning("Input tree branch length ignored");
    _data.is_length_set=false;
  }

  if(is_leaf()) return;
  child1()->clear_length(is_warn);
  child2()->clear_length(is_warn);
}



void
bi_tree_node::
make_index(std::vector<bi_tree_node*>& node_index,
           std::vector<const bi_tree_node*>& leaf_node_index){

  if(label().empty()){
    std::ostringstream oss;
    oss << "node" << node_index.size();
    _data.label=oss.str();
  }

  this->_data.node_id=node_index.size();
  node_index.push_back(this);

  if(is_leaf()) {
    this->_data.leaf_id=leaf_node_index.size();
    leaf_node_index.push_back(this);
  } else {
    child1()->make_index(node_index,leaf_node_index);
    child2()->make_index(node_index,leaf_node_index);
  }
}



bi_tree::
bi_tree(const char* tree_string)
  : _root(new bi_tree_node(newick_tree_parser(tree_string).root())) {
  ctor_setup_tree_index();
}



bi_tree::
bi_tree(const bi_tree& s)
  : _root(new bi_tree_node(s.root())) {
  ctor_setup_tree_index();
}



static
std::ostream& operator<<(std::ostream& os,
                         const bi_tree_node* b){

  if(b) {
    if(! b->is_leaf()){
      os << "(" << b->child1() << "," << b->child2() << ")";
    }
    os << b->label();
    if(b->is_length_set()) os << ":" << b->length();
  }
  return os;
}



std::ostream& operator<<(std::ostream& os,
                         const bi_tree& t){
  return os << t.root() << ";";
}



void
bi_tree::
ctor_setup_tree_index(){
  root()->make_index(_node_index,_leaf_node_index);

  for(unsigned i(0);i<node_size();++i){
    _node_label_map[node(i)->label()]=node(i);
  }
}



void
bi_tree::
clear_branch_lengths(const bool is_warn){
  root()->clear_length(is_warn);
}

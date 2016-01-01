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

// $Id: newick_tree_parser.cc 1207 2008-05-14 19:22:13Z ctsa $

/// \file

#include "newick_tree_parser.h"

#include <cstdlib>

#include <ostream>


namespace NEWICK_PARSER {

  static
  void
  parse_error(const char* err,
              const char* const s_start,
              const char* s){
    static const std::string h(std::string("newick tree parser: ")+err+
                               std::string("\n...original tree string:\n")+s_start+
                               std::string("\n...parse error after:\n")+std::string(s_start,s-s_start));
    throw substk_exception(h.c_str());
  }



  enum {
    START='(',
    STOP=')',
    JOIN=',',
    LABELD=':',
    END=';'
  };



  inline
  static
  bool
  is_symbol(const char c){
    switch(c){
    case START:
    case STOP:
    case JOIN:
    case LABELD:
    case END:    return true;
    default:     return false;
    }
  }



  inline
  static
  const char*
  next_token(const char*& s){
    while(isspace(*s)) ++s;
    return s;
  }



  static
  void
  parse_label(const char* const s_start,
              const char*& s,
              std::string& label){
    if(! label.empty()) parse_error("malformed tree string",s_start,s);

    const char* start(next_token(s));
    const char* stop(start);
    for(;;++s){
      if       (*s!=' ' && (is_symbol(*s) || isspace(*s) || *s == 0)){
        label=std::string(start,stop-start);
        return;
      } else if(*s!=' '){
        stop=s+1;
      }
    }
  }

  static
  void
  parse_value(const char* const s_start,
              const char*& s,
              newick_tree_node* b){
    char* s_end(0);
    b->value=strtod(next_token(s),&s_end);
    if(s==s_end) parse_error("can't parse node value to floating point",s_start,s);
    b->is_value_set=true;
    s=s_end;
  }

  static
  void
  parse_label_value(const char* const s_start,
                    const char*& s,
                    newick_tree_node* b){
    parse_label(s_start,s,b->label);
    if(*next_token(s)==LABELD){ parse_value(s_start,++s,b); }
  }



  static
  newick_tree_node*
  parse_subtree(const char* const s_start,
                const char*& s){

    newick_tree_node* n(new newick_tree_node());

    if(*next_token(s)==START){
      do{
        n->add_child(parse_subtree(s_start,++s));
      } while(*next_token(s)==JOIN);

      if(*s!=STOP) parse_error("expected ')' to close paren block",s_start,s);
      ++s;
    }
    parse_label_value(s_start,s,n);
    return n;
  }



  static
  newick_tree_node*
  parse(const char* s){
    const char* const s_start(s);
    newick_tree_node* n(parse_subtree(s_start,s));
    if(*next_token(s)==0) parse_error("unexpected termination (possibly no closing semicolon)",s_start,s);
    if(*s!=END)           parse_error("malformed tree string",s_start,s);
    return n;
  }
}



newick_tree_parser::
newick_tree_parser(const char* tree_string)
  : _root(NEWICK_PARSER::parse(tree_string)) {}



static
std::ostream& operator<<(std::ostream& os,
                         const newick_tree_node* b){

  if(b==0) return os;

  const unsigned nc(b->child_size());
  if(nc){
    os << "(";
    for(unsigned i(0);i<nc;++i){
      if(i) os << ",";
      os << b->child(i);
    }
    os << ")";
  }
  os << b->label;
  if(b->is_value_set) os << ":" << b->value;
  return os;
}



std::ostream& operator<<(std::ostream& os,
                         const newick_tree_parser& t){
  return os << t.root() << ";";
}


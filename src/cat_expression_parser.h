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

// $Id: newick_tree_parser.h 931 2007-10-22 19:07:34Z ctsa $

/// \file

#ifndef __CAT_EXPRESSION_PARSER_H
#define __CAT_EXPRESSION_PARSER_H


#include "util/general/uncopyable.h"

#include <string>
#include <vector>



///
struct cat_expression_parser : private uncopyable {

  typedef cat_expression_parser self_t;

  cat_expression_parser(const char* ce);

  typedef std::pair<std::string,std::string> param_set_mapping_rule;
  typedef std::vector<param_set_mapping_rule> param_set_mappings;
  typedef std::pair<std::string,param_set_mappings> param_set_definition;
  typedef std::vector<param_set_definition> param_set_definitions;
  typedef std::pair<std::string,param_set_definitions> seq_cat_definition;
  typedef std::vector<seq_cat_definition> seq_cat_definitions;

  typedef std::vector<std::string> data_set_seq_cats;
  typedef std::vector<std::string> data_class_labels;
  typedef std::pair<std::string,data_class_labels> data_set_definition;
  typedef std::pair<data_set_definition,data_set_seq_cats> data_set_cat_map_rule;
  typedef std::vector<data_set_cat_map_rule> data_set_cat_mapping;

  typedef std::pair<seq_cat_definitions,data_set_cat_mapping> cat_expression;

  const cat_expression& ce() const { return _ce; }

private:
  cat_expression _ce;
};

#endif

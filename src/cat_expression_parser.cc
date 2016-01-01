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


#include "cat_expression_parser.h"
#include "substk_exception.h"



/// bogus class is here to make improved parse error messages easy by
/// keeping s_in & s handy
///
struct ce_parser_internal {

  typedef cat_expression_parser cep;


  ce_parser_internal(const char* s_init,
                     cep::cat_expression& ce)
    : s_in(s_init), s(s_init){

    parse_seq_cat_definitions(ce.first);
    parse_data_set_cat_mapping(ce.second);
  }


private:


  void
  parse_error(const char* err){
    static const std::string h(std::string("cat expression parser: ")+err+
                               std::string("\n...original cat expression string:\n")+s_in+
                               std::string("\n...parse error after:\n")+std::string(s_in,s-s_in));
    throw substk_exception(h.c_str());
  }


  const char*
  next_token(){
    while(isspace(*s)) ++s;
    return s;
  }


  enum ce_syntax {
    param_set_label_split=':',
    param_set_mapping_rule_delimit=',',
    param_set_mappings_open='{',
    param_set_mappings_close='}',
    param_set_definitions_open='{',
    param_set_definitions_close='}',
    data_class_label_delimit=',',
    data_class_labels_open='{',
    data_class_labels_close='}',
    data_set_cat_map_rule_open='[',
    data_set_cat_map_rule_close=']'
  };


  static
  bool
  is_symbol(const char c){
    switch(c){
    case param_set_label_split:
    case param_set_mapping_rule_delimit:
    case param_set_definitions_open:
    case param_set_definitions_close:
    case data_set_cat_map_rule_open:
    case data_set_cat_map_rule_close: return true;
    default:                          return false;
    }
  }


  void
  parse_label(std::string& label){
    if(! label.empty()) parse_error("malformed cat expression string");

    const char* start(next_token());
    while(! (is_symbol(*s) || isspace(*s) || *s == 0)) ++s;
    label=std::string(start,s-start);

    if(label.empty()) parse_error("empty label values not allowed");
  }


  void
  parse_param_set_mapping_rule(cep::param_set_mapping_rule& psmr){
    parse_label(psmr.first);
    if(*next_token()==param_set_label_split) {
      ++s;
      parse_label(psmr.second);
    }
  }


  void
  parse_param_set_mappings(cep::param_set_mappings& psm){
    if(*next_token()!=param_set_mappings_open) parse_error("expected opening brace of param_set_mappings");
    do {
      ++s;
      psm.resize(psm.size()+1);
      parse_param_set_mapping_rule(psm.back());
    } while(*next_token()==param_set_mapping_rule_delimit);
    if(*s!=param_set_mappings_close) parse_error("expected closing brace of param_set_mappings");
    ++s;
  }


  void
  parse_param_set_definition(cep::param_set_definition& psd){
    parse_label(psd.first);
    if(psd.first.size() != 3) parse_error((std::string("invalid param_set_code: ")+psd.first).c_str());
    parse_param_set_mappings(psd.second);
  }


  void
  parse_param_set_definitions(cep::param_set_definitions& psds){
    if(*next_token()!=param_set_definitions_open) parse_error("expected opening brace of param set definitions");
    ++s;
    do{
      psds.resize(psds.size()+1);
      parse_param_set_definition(psds.back());
      if(next_token()==0) parse_error("can't find closing brace of param set definitions");
    } while(*s!=param_set_definitions_close);
    ++s;
  }


  void
  parse_seq_cat_definition(cep::seq_cat_definition& scd){
    parse_label(scd.first);
    parse_param_set_definitions(scd.second);
  }


  void
  parse_seq_cat_definitions(cep::seq_cat_definitions& scd){
    while(true){
      next_token();
      if(*s==data_set_cat_map_rule_open || *s==0) break;
      scd.resize(scd.size()+1);
      parse_seq_cat_definition(scd.back());
    }
  }


  void
  parse_data_set_seq_cats(cep::data_set_seq_cats& dssc){
    while(*next_token()!=data_set_cat_map_rule_close) {
      if(is_symbol(*s) || *s==0){
        parse_error("expected closing brace of data set cat map rule");
      }
      dssc.resize(dssc.size()+1);
      parse_label(dssc.back());
    }
  }


  void
  parse_data_class_labels(cep::data_class_labels& dcl){
    if(*next_token()!=data_class_labels_open) {
      parse_error("expected opening brace of data class labels");
    }
    do {
      ++s;
      dcl.resize(dcl.size()+1);
      parse_label(dcl.back());
    } while(*next_token()==data_class_label_delimit);
    if(*s!=data_class_labels_close) {
      parse_error("expected closing brace of data class labels");
    }
    ++s;
  }


  void
  parse_data_set_definition(cep::data_set_definition& dsd){
    parse_label(dsd.first);
    parse_data_class_labels(dsd.second);
  }


  void
  parse_data_set_cat_map_rule(cep::data_set_cat_map_rule& dsmr){
    if(*next_token()!=data_set_cat_map_rule_open){
      parse_error("expected opening brace of data_set_cat_map_rule");
    }
    ++s;

    parse_data_set_definition(dsmr.first);
    parse_data_set_seq_cats(dsmr.second);

    if(*next_token()!=data_set_cat_map_rule_close) {
      parse_error("expected closing brace of data_set_cat_map_rule");
    }
    ++s;
  }


  void
  parse_data_set_cat_mapping(cep::data_set_cat_mapping& dsm){
    while(*next_token() != 0){
      dsm.resize(dsm.size()+1);
      parse_data_set_cat_map_rule(dsm.back());
    }
  }

  const char* const s_in;
  const char* s;
};



cat_expression_parser::
cat_expression_parser(const char* ce_string){
  ce_parser_internal(ce_string,_ce);
}

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

// $Id: cat_manager.h 1222 2008-05-22 23:10:06Z ctsa $

/// \file

#ifndef __CAT_MANAGER_H
#define __CAT_MANAGER_H

#include "cat_expression_parser.h" /// \todo include here only for private method arg, get rid of somehow...
#include "cat_info.h"
#include "cat_type.h"
#include "name_id_lup.h"
#include "param_composite.h"
#include "param_init_type.h"
#include "prob_gtor_param.h"
#include "simple_util.h"
#include "subs_ml_model_comex.h"
#include "subs_ml_ptol.h"
#include "subs_ml_types.h"

#include <algorithm>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>



// non-persistent initialization options
struct cat_manager_options {
  cat_manager_options()
    : is_lock_cat_prob(true),
      cat_model_str() {}

  bool is_lock_cat_prob;
  std::string cat_model_str;
};



typedef std::map<std::string,unsigned> data_class_label_assignment_map_type;



struct cat_manager : public param_composite {

  typedef param_composite base_t;
  typedef cat_manager self_t;

  cat_manager(const cat_manager_options& cmo,
              const cat_manager_sml_share& smls);

  cat_manager(std::istream& is,
              const cat_manager_sml_share& smls);

  cat_manager(const self_t& s,
              const cat_manager_sml_share& smls) :
    base_t(s),
    _smls(smls),
    _ptor(s.ptor().clone()),
    _data(s._data) { register_param(ptor()); }

private:
  self_t& operator=(const self_t&);
public:

  unsigned
  typed_cat_no_from_cat_no(const unsigned cat_no,
                           const CAT_PARAM_TYPE::index_t pt,
                           const CAT_MIX_TYPE::index_t mt = CAT_MIX_TYPE::EITHER) const {
    /// is there a way to enforce this limitation at compile time??
    /// could we pull off templating this method?
    /// \todo make this a template on pt
    if(is_branch_cat_param_type(pt)){
      die("disallowed call to typed_cat_no_from_cat_no");
    }

    static const unsigned branch_id(0);
    return typed_cat_no_from_cat_no_y_branch_id(cat_no,branch_id,pt,mt);
  }

  unsigned
  typed_cat_no_from_cat_no_y_branch_cat_set(const unsigned cat_no,
                                            const unsigned branch_cat_set,
                                            const CAT_PARAM_TYPE::index_t pt,
                                            const CAT_MIX_TYPE::index_t mt = CAT_MIX_TYPE::EITHER) const {
    const unsigned branch_id(_data.branch_cat_set_first_branch_id[cat_no][branch_cat_set]);
    return  typed_cat_no_from_cat_no_y_branch_id(cat_no,branch_id,pt,mt);
  }

  unsigned
  typed_cat_no_from_cat_no_y_branch_id(const unsigned cat_no,
                                       const unsigned branch_id,
                                       const CAT_PARAM_TYPE::index_t pt,
                                       const CAT_MIX_TYPE::index_t mt = CAT_MIX_TYPE::EITHER) const {
    if(mt == CAT_MIX_TYPE::EITHER){
      if(CAT_PARAM_TYPE::is_mixable[pt])  die("invalid call to typed_cat_size");
      return std::max(typed_cat_no_from_cat_no_y_branch_id(cat_no,branch_id,pt,CAT_MIX_TYPE::SITE),
                      typed_cat_no_from_cat_no_y_branch_id(cat_no,branch_id,pt,CAT_MIX_TYPE::GROUP));
    } else {
      return _data.typed_cat_no[cat_no][branch_id][typed_cat_index(pt,mt)];
    }
  }

#if 0
  unsigned
  typed_cat_no_from_branch_typed_cat_no_y_branch_id(const unsigned branch_typed_cat_no,
                                                    const unsigned,
                                                    const CAT_PARAM_TYPE::index_t,
                                                    const CAT_MIX_TYPE::index_t mt = CAT_MIX_TYPE::EITHER) const {
    /// finish impl
    return branch_typed_cat_no;
  }
#endif

  unsigned
  typed_cat_size(const CAT_PARAM_TYPE::index_t pt,
                 const CAT_MIX_TYPE::index_t mt = CAT_MIX_TYPE::EITHER) const {
    if(mt == CAT_MIX_TYPE::EITHER){
      if(CAT_PARAM_TYPE::is_mixable[pt])  die("invalid call to typed_cat_size");
      return std::max(typed_cat_size(pt,CAT_MIX_TYPE::SITE),typed_cat_size(pt,CAT_MIX_TYPE::GROUP));
    } else {
      return _data.typed_cat_size[typed_cat_index(pt,mt)];
    }
  }

  //
  unsigned
  group_cat_size() const {
    return _data.cat_to_group.back()+1;
  }

  unsigned
  group_cat_no(const unsigned cat_no) const {
    return _data.cat_to_group[cat_no];
  }

  unsigned
  cat_size() const {
    return _data.typed_cat_no.dim1();
  }

  unsigned
  assigned_data_set_size() const {
    /// \todo ads_fix
    return 1;
#if 0
    const data_class_label_assignment_map_type& dc_ads_map(data_class_label_assignment_map());
    if(dc_ads_map.empty()) {
      return 1;
    } else {
      // more stuff....
    }
#endif
  }

  /// \brief true for each cat_no being used in the mixture tied to
  /// adset,
  ///
  /// note that this refers to the parameter-independent relationship,
  /// so the adset->cat relationship is not set to false just because
  /// p=0 in a particlular model parameterization
  ///
  /// \todo ads_fix
  void
  is_adset_using_cat(bool* iac,
                     const unsigned /*adset_id*/) const {

    const unsigned n_cats(cat_size());
    for(unsigned c(0);c<n_cats;++c){
      iac[c]=true;
    }
  }

  const data_class_label_assignment_map_type&
  data_class_label_assignment_map() const { return _data.data_class_label_assignment_map; }

  void store_state(std::ostream& os) const;

  template <typename RandomAccessIterator>
  void
  cat_pdistro(RandomAccessIterator x) const {
    const unsigned cs(cat_size());
    const prob_t* cp(ptor().pdistro());
    for(unsigned i(0);i<cs;++i){ *(x+i) = cp[i]; }
    pdistro_norm(x,x+cs);
#ifdef DEBUG
    pdistro_check(x,cs,SUBS_ML_PTOL);
#endif
  }

  /// \todo ads_fix
  template <typename RandomAccessIterator>
  void
  adset_cat_pdistro(RandomAccessIterator x,
                    const unsigned /*adset_id*/) const {
    cat_pdistro(x);
  }


  template <typename RandomAccessIterator>
  void
  group_cat_pdistro(RandomAccessIterator x) const {
    const unsigned cs(cat_size());
    const unsigned gcs(group_cat_size());
    std::fill(x,x+gcs,0.);

    const prob_t* cp(ptor().pdistro());
    for(unsigned i(0);i<cs;++i){
      const unsigned gc(group_cat_no(i));
      *(x+gc) += cp[i];
    }
    pdistro_norm(x,x+gcs);
#ifdef DEBUG
    pdistro_check(x,gcs,SUBS_ML_PTOL);
#endif
  }

  /// \todo ads_fix
  template <typename RandomAccessIterator>
  void
  adset_group_cat_pdistro(RandomAccessIterator x,
                          const unsigned /*adset_id*/) const {
    group_cat_pdistro(x);
  }

  /// \todo pdistro's should be reported as irv_t's too, via a template version
  /// of this function:
  ///
  void
  typed_cat_pdistro(prob_t* tcp,
                    const CAT_PARAM_TYPE::index_t pt,
                    const CAT_MIX_TYPE::index_t mt = CAT_MIX_TYPE::EITHER) const {
    if(is_branch_cat_param_type(pt)){
      die("disallowed call to typed_cat_no_from_cat_no");
    }

    static const unsigned branch_id(0);
    branch_typed_cat_pdistro(tcp,branch_id,pt,mt);
  }

  void
  branch_typed_cat_pdistro(prob_t* tcp,
                           const unsigned branch_id,
                           const CAT_PARAM_TYPE::index_t pt,
                           const CAT_MIX_TYPE::index_t mt = CAT_MIX_TYPE::EITHER) const;

  const std::string&
  cat_label(const unsigned cat_no) const;

  const std::string&
  typed_cat_label(const unsigned typed_cat_no,
                  const CAT_PARAM_TYPE::index_t pt,
                  const CAT_MIX_TYPE::index_t mt = CAT_MIX_TYPE::EITHER) const;

  void
  reset_param(const PARAM_INIT_TYPE::index_t pinit);

  void
  set_cat_pdistro(const smlfloat* x){
    set_param_state(x);
  }

  unsigned
  mut_model_norm_set_size() const {
    return 1+*std::max_element(_data.mut_model_set_by_mut_model_cat.begin(),
                               _data.mut_model_set_by_mut_model_cat.end());
  }

  std::vector<unsigned>::const_iterator
  mut_model_norm_set_from_mut_model_cat() const {
    return _data.mut_model_set_by_mut_model_cat.begin();
  }

  unsigned
  mut_rate_norm_set_size(const CAT_MIX_TYPE::index_t mt) const {
    return 1+*std::max_element(_data.mut_rate_set_by_mut_rate_cat[mt].begin(),
                               _data.mut_rate_set_by_mut_rate_cat[mt].end());
  }

  std::vector<unsigned>::const_iterator
  mut_rate_norm_set_from_mut_rate_cat(const CAT_MIX_TYPE::index_t mt) const {
    return _data.mut_rate_set_by_mut_rate_cat[mt].begin();
  }

  unsigned
  sel_strength_norm_set_size(const CAT_MIX_TYPE::index_t mt) const {
    return 1+*std::max_element(_data.sel_strength_set_by_sel_strength_cat[mt].begin(),
                               _data.sel_strength_set_by_sel_strength_cat[mt].end());
  }

  std::vector<unsigned>::const_iterator
  sel_strength_norm_set_from_sel_strength_cat(const CAT_MIX_TYPE::index_t mt) const {
    return _data.sel_strength_set_by_sel_strength_cat[mt].begin();
  }

#if 0
  unsigned
  mut_model_norm_set_from_time_cat(const unsigned time_cat_no) const {
    assert(time_cat_no<typed_cat_size(CAT_PARAM_TYPE::TIME));
    return _data.mut_model_set_by_time_cat[time_cat_no];
  }
#endif

  //  bool is_missing_obs_partition() const;

  unsigned branch_size() const;

  const std::string& branch_label(const unsigned branch_id) const;

  void branch_cat_set_label(const unsigned cat_no,
                            const unsigned branch_cat_set,
                            std::string& label) const;

  unsigned branch_cat_set_size(const unsigned cat_no) const {
    return _data.branch_cat_set_first_branch_id[cat_no].size();
  }

  /// \todo branch_cat_set should be identifiable from branch_id alone,
  /// simplify this or add better doc of what branch_cat_set means
  ///
  unsigned get_branch_cat_set(const unsigned cat_no,
                              const unsigned branch_id) const {
    return _data.branch_cat_set_no[cat_no][branch_id];
  }

  bool is_missing_obs_partition() const;


  void report(std::ostream& os) const;

private:

  unsigned&
  typed_cat_no_from_cat_no_y_branch_id(const unsigned cat_no,
                                       const unsigned branch_id,
                                       const CAT_PARAM_TYPE::index_t pt,
                                       const CAT_MIX_TYPE::index_t mt) {
    return _data.typed_cat_no[cat_no][branch_id][typed_cat_index(pt,mt)];
 }

  void init_cat_info();

  void
  add_param_set_mapping_rule(const CAT_PARAM_TYPE::index_t& pt,
                             const CAT_MIX_TYPE::index_t& mt,
                             const std::string& label,
                             unsigned* typed_cat_no,
                             bool* is_typed_cat_no);

  void
  process_seq_cat_definitions(const cat_expression_parser::seq_cat_definitions& scd);

  void
  parse_cat_expression(const std::string&);

  prob_gtor_param& ptor() { return *_ptor; }
  const prob_gtor_param& ptor() const { return *_ptor; }

  bool
  is_branch_cat_param_type(const CAT_PARAM_TYPE::index_t pt) const {
    using namespace CAT_PARAM_TYPE;
    switch(pt){
    case MUT_RATE:
    case MUT_MODEL:
    case SEL_STRENGTH:
    case SEL_MATRIX:    return true;
    default:            return false;
    }
  }

  struct auto_copy {
    auto_copy() :
      typed_cat_size(CAT_TYPE_SIZE,0),
      typed_cat_label_lup(CAT_TYPE_SIZE),
      mut_rate_set_by_mut_rate_cat(CAT_MIX_TYPE::SIZE),
      sel_strength_set_by_sel_strength_cat(CAT_MIX_TYPE::SIZE) {}

    std::vector<unsigned> cat_to_group;
    data_class_label_assignment_map_type data_class_label_assignment_map;
    simple_init_matrix3d<unsigned> typed_cat_no;
    std::vector<unsigned> typed_cat_size;
    name_id_lup cat_label_lup;
    std::vector<name_id_lup> typed_cat_label_lup;
    std::vector<std::vector<unsigned> > branch_cat_set_first_branch_id;
    simple_init_matrix<unsigned> branch_cat_set_no;
    std::vector<unsigned> mut_model_set_by_mut_model_cat;
    std::vector<std::vector<unsigned> > mut_rate_set_by_mut_rate_cat;
    std::vector<std::vector<unsigned> > sel_strength_set_by_sel_strength_cat;
    //    std::vector<unsigned> mut_model_set_by_time_cat;
  };


  const cat_manager_sml_share _smls;
  std::auto_ptr<prob_gtor_param> _ptor;
  auto_copy _data;
};


#endif

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

// $Id: subs_ml_model.h 1148 2008-03-11 00:49:49Z ctsa $

/// \file

#ifndef __SUBS_ML_MODEL_H
#define __SUBS_ML_MODEL_H


#include "audit_info.h"
#include "param_composite.h"
#include "param_init_type.h"
#include "rate_gtor_options.h"
#include "simple_util.h"
#include "subs_ml_model_comex.h"

#include <iosfwd>
#include <memory>


struct bg_gtor;
struct bi_tree;
struct cat_manager;
struct notifier;
struct obs_info;
struct param_composite;
struct rate_gtor;
struct root_gtor;
struct site_data;
struct site_data_fastlup;
struct subs_ml_model_init_options;
struct subs_ml_model_min_options;
struct time_gtor;



/// \brief parameters and behavior of ml model
///
struct subs_ml_model : public param_composite {

  typedef param_composite base_t;
  typedef subs_ml_model self_t;

  explicit
  subs_ml_model(const subs_ml_model_init_options& sml_opt);

  subs_ml_model(std::istream& is,
                const audit_info& ai);

  subs_ml_model(const self_t& s);

  ~subs_ml_model();

private:
  self_t& operator=(const self_t&);
public:

  RATE_GTOR_MODEL::index_t rate_gtor_model() const { return _data.ragm; }

  const time_gtor& get_time_gtor() const { return *_time_gtor; }
  const rate_gtor& get_rate_gtor() const { return *_rate_gtor; }
  const bg_gtor& get_bg_gtor() const { return *_bg_gtor; }
  const root_gtor& get_root_gtor() const { return *_root_gtor; }
  const cat_manager& get_cat_manager() const { return *_cm; }
  const bi_tree& tree() const { return *_tree; }

  subs_ml_model_min_options& opt() {return *_opt;}
  const subs_ml_model_min_options& opt() const {return *_opt;}

  const notifier& bg_gtor_notifier() const;
  const notifier& obs_info_notifier() const;

  unsigned state_size() const;

  unsigned submodel_size() const;

  unsigned submodel_state_size(unsigned s) const;

  void store_state(std::ostream& os) const;

  smlfloat param_penalty() const;

  /// \brief report information processed from model parameters
  void report(std::ostream& os) const;

  /// \brief time subsection of report
  void report_time(std::ostream& os) const;

  const audit_info& get_audit_info() const { return _data.ai; }

  const prob_t * const * obs_seq_state_distro_cat(const unsigned cat_no) const;

  const prob_t * const * obs_seq_state_distro_obs_cat(const unsigned obs_cat_no) const;

  const prob_t * obs_state_distro_cat(const unsigned cat_no) const;

  const prob_t * obs_state_distro_obs_cat(const unsigned obs_cat_no) const;

  /// \brief trained parameters initialized from data -- should be
  /// called for new models only
  ///
  /// 1) set the root distro to that of the present day data
  /// 2) set assigned cat priors
  /// 3) match assigned cat type to data cats
  ///
  void init_trained_data_dependencies(const site_data& sd);

  /// \brief untrained parameters dependent on data -- should be called
  /// to initialize any model
  ///
  /// 1) set the background distribution to that of the present day data
  ///
  void init_untrained_data_dependencies(const site_data& sd);

  void norm();

  /// returns is_delta_lnL
  ///
  bool
  update_model_cat_expectation(const site_data_fastlup& sdf,
                               const smlfloat ntol);

  void update_model_obs_distro(const site_data_fastlup& sdf);

  void check_cat_post_prob_data_match(const site_data_fastlup& sdf);

  // very temporary method to access data
  //
  bool is_gcp() const {
    return _data.group_cat_post_prob.dim1() != 0;
  }

  const prob_t* const * get_gcp() const {
    return _data.group_cat_post_prob.ptr();
  }

#if 0
  std::vector<std::string>::const_iterator get_glabel() const {
    return static_cast<const std::vector<std::string>&>(_data.group_label).begin();
  }
#endif

  bool& is_cat_mstep() { return _data.is_cat_mstep; }
  bool is_cat_mstep() const { return _data.is_cat_mstep; }

  /// \todo convert this to friend access...
  root_gtor& get_root_gtor_nonconst(){ return *_root_gtor; }

private:
  time_gtor& get_time_gtor_nonconst(){ return *_time_gtor; }
  rate_gtor& get_rate_gtor_nonconst(){ return *_rate_gtor; }
  bg_gtor& get_bg_gtor_nonconst() { return *_bg_gtor; }
  cat_manager& get_cat_manager_nonconst() { return *_cm; }

  obs_info& obs_nonconst() { return *_obs; }
  const obs_info& obs() const { return *_obs; }

  bi_tree& tree_nonconst() { return *_tree; }

  friend bi_tree& time_gtor_sml_share::tree_nonconst();

  void set_rate_gtor_model(const subs_ml_model_init_options& sml_opt);

  void reset_param(const PARAM_INIT_TYPE::index_t pinit);

  void check_param() const;

  void register_param();

  // notify components which read info from sml model (through
  // subs_ml_model_comex) that the parameter state has been updated:
  //
  virtual void post_write_param_state();

  void set_obs_cat_seq_state_counts(smlfloat const * const * const * c);
  void set_obs_cat_seq_state_distro(prob_t const * const * const * c);

  void clear_post_prob_data();

  //////
  struct auto_copy {
    auto_copy(const audit_info& _ai,
              const RATE_GTOR_MODEL::index_t _ragm = RATE_GTOR_MODEL::NONE)
      : ai(_ai), ragm(_ragm), is_cat_mstep(false) {}

    const audit_info ai;
    RATE_GTOR_MODEL::index_t ragm;

    simple_init_matrix<prob_t> group_cat_post_prob;
    std::vector<std::string> group_label;
    bool is_cat_mstep;
  };

  auto_copy _data;
  std::auto_ptr<bi_tree> _tree;
  std::auto_ptr<cat_manager> _cm;
  std::auto_ptr<time_gtor> _time_gtor;
  std::auto_ptr<obs_info> _obs;
  std::auto_ptr<bg_gtor> _bg_gtor;
  std::auto_ptr<rate_gtor> _rate_gtor;
  std::auto_ptr<root_gtor> _root_gtor;
  std::auto_ptr<subs_ml_model_min_options> _opt;
};

#endif

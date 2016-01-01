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

// $Id: subs_ml_model.cc 1210 2008-05-20 23:17:36Z ctsa $

/// \file

#include "bg_gtor.h"
#include "rate_gtor.h"
#include "rate_gtor_factory.h"
#include "root_gtor.h"
#include "report_util.h"
#include "subs_ml_model.h"
#include "subs_ml_model_init_options.h"
#include "subs_ml_model_min_options.h"
#include "util/general/die.h"
#include "util/general/io_util.h"
#include "util/general/log.h"
#include "util/general/metatags.h"
#include "util/math/prob_util_io.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>


const unsigned MODEL_FILE_VERSION(10);


const char* const section_id = "subs_ml_model";
const char* const end_label = "END";

const char* const gcp_label = "group_cat_prob";



void
subs_ml_model::
store_state(std::ostream& os) const {

  os << section_id << " model_file_version " << MODEL_FILE_VERSION << "\n";
  get_audit_info().store_state(os);
  os << section_id << " tree " << tree() << "\n";
  get_cat_manager().store_state(os);
  get_time_gtor().store_state(os);
  os << section_id << " rate_gtor_model " << rate_gtor_model() << "\n";
  obs().store_state(os);
  get_bg_gtor().store_state(os);
  get_rate_gtor().store_state(os);
  get_root_gtor().store_state(os);
  opt().store_state(os);

  if(is_gcp()){
    assert(_data.group_label.size() == _data.group_cat_post_prob.dim1());
    const unsigned n_cats(get_cat_manager().cat_size());
    const unsigned n_groups(_data.group_label.size());
    os << gcp_label << " n_groups " << n_groups << "\n";

    for(unsigned i(0);i<n_groups;++i){
      os << gcp_label << " " << _data.group_label[i];
      for(unsigned c(0);c<n_cats;++c){
        os << " " << _data.group_cat_post_prob[i][c];
      }
      os << "\n";
    }
    os << gcp_label << " " << end_label << "\n";
  }

  os << section_id << " " << end_label << "\n";
}



static
void
load_die(std::istream& is,
         const char* const sid,
         const std::string& tmpstr)  NORETURN_TAG;

static
void
load_die(std::istream& is,
         const char* const sid,
         const std::string& tmpstr){

  log_os << "ERROR:: -- Invalid model file format. --\n"
         << "ERROR:: expected_section: " << sid << "\n"
         << "ERROR:: section: " << tmpstr << "\n"
         << "\n";
  log_os << is.rdbuf();
  abort();
}



static
void
load_advance(std::istream& is,
             const char* const second = 0){

  std::string tmpstr;
  is >> tmpstr;
  if(tmpstr != section_id) load_die(is,section_id,tmpstr);
  is >> tmpstr;
  if(second && tmpstr != second) load_die(is,second,tmpstr);
}

static
void
load_advance_gcp(std::istream& is,
                 const char* const second = 0){

  std::string tmpstr;
  is >> tmpstr;
  if(tmpstr != gcp_label) load_die(is,gcp_label,tmpstr);
  if(second){
    is >> tmpstr;
    if(tmpstr != second) load_die(is,second,tmpstr);
  }
}


subs_ml_model::
subs_ml_model(const subs_ml_model_init_options& sml_opt)
  : base_t(), _data(sml_opt.ai,sml_opt.ragm),
    _tree(new bi_tree(sml_opt.tree_buffer.c_str())),
    _cm(new cat_manager(sml_opt.catman_opt,cat_manager_sml_share(*this))),
    _time_gtor(new time_gtor(sml_opt.tgo,time_gtor_sml_share(*this))),
    _opt(new subs_ml_model_min_options()){

  assert(rate_gtor_model() != RATE_GTOR_MODEL::NONE);

  {
    const unsigned n_obs_cats(get_cat_manager().typed_cat_size(CAT_PARAM_TYPE::OBS));
    const unsigned n_seqs(tree().leaf_size());
    const unsigned n_states(state_size());
    _obs.reset(new obs_info(n_obs_cats,n_seqs,n_states));
  }

  _bg_gtor.reset(new bg_gtor(sml_opt.bopt,bg_gtor_sml_share(*this)));
  set_rate_gtor_model(sml_opt);
  _root_gtor.reset(new root_gtor(sml_opt.rogm,root_gtor_sml_share(*this)));
  register_param();

  reset_param(sml_opt.pinit);
}



subs_ml_model::
subs_ml_model(const self_t& s)
  :  base_t(s), _data(s._data),
     _tree(new bi_tree(s.tree())),
     _cm(new cat_manager(s.get_cat_manager(),cat_manager_sml_share(*this))),
     _time_gtor(new time_gtor(s.get_time_gtor(),time_gtor_sml_share(*this))),
     _obs(new obs_info(s.obs())),
     _bg_gtor(new bg_gtor(s.get_bg_gtor(),bg_gtor_sml_share(*this))),
     _rate_gtor(s.get_rate_gtor().clone(rate_gtor_sml_share(*this))),
     _root_gtor(new root_gtor(s.get_root_gtor(),root_gtor_sml_share(*this))),
     _opt(new subs_ml_model_min_options(s.opt()))
{ register_param(); }



/// \brief 'unserializes' subs_ml_model state from an input stream
///
/// note that the sml model file is actually laid out in the order of
/// initialization dependency, and in an unacceptably brittle way -- So:
///
/// \todo separate stream parsing code from sml restoration
///
/// \todo set up checks on initialization order, so that mistakes
/// don't occur as new dependencies are introduced among sml
/// components
///
subs_ml_model::
subs_ml_model(std::istream& is,
              const audit_info& ai) : base_t(), _data(ai) {

  std::string dummy;

  {  // check file version
    unsigned tmpu;
    load_advance(is,"model_file_version");
    is >> tmpu;
    if(tmpu != MODEL_FILE_VERSION){
      log_os << "ERROR:: bad model file version number. expected: " << MODEL_FILE_VERSION << "\n";
      exit(EXIT_FAILURE);
    }
  }

  // throw away old audit_info
  get_audit_info().skip_state(is);

  // get the phylogenetic tree
  load_advance(is,"tree");
  is >> dummy;
  _tree.reset(new bi_tree(dummy.c_str()));
  tree_nonconst().clear_branch_lengths();

  _cm.reset(new cat_manager(is,cat_manager_sml_share(*this)));
  _time_gtor.reset(new time_gtor(is,time_gtor_sml_share(*this)));

  // ragm is needed by all other components, so it comes first in the file:
  //
  load_advance(is,"rate_gtor_model");
  _data.ragm=static_cast<RATE_GTOR_MODEL::index_t>(read_int(is));
  assert(rate_gtor_model() != RATE_GTOR_MODEL::NONE);

  const unsigned n_obs_cats(get_cat_manager().typed_cat_size(CAT_PARAM_TYPE::OBS));
  const unsigned n_seqs(tree().leaf_size());
  const unsigned n_states(state_size());

  _obs.reset(new obs_info(n_obs_cats,n_seqs,n_states,is));

  _bg_gtor.reset(new bg_gtor(is,bg_gtor_sml_share(*this)));

  const subs_ml_model_init_options sml_opt(_data.ragm,ai);
  set_rate_gtor_model(sml_opt);
  get_rate_gtor_nonconst().load_state(is);

  _root_gtor.reset(new root_gtor(is,root_gtor_sml_share(*this)));

  _opt.reset(new subs_ml_model_min_options(is));

  /// \todo put this gcp business into a separate object:
  ///
  is >> dummy;
  if(dummy==gcp_label){
    const unsigned n_cats(get_cat_manager().cat_size());

    is >> dummy;
    static const char * const ng_label = "n_groups";
    if(dummy != ng_label){
      load_die(is,ng_label,dummy);
    }

    unsigned n_groups;
    is >> n_groups;

    _data.group_label.resize(n_groups);
    _data.group_cat_post_prob.init(n_groups,n_cats);

    for(unsigned i(0);i<n_groups;++i){
      load_advance_gcp(is);
      is >> _data.group_label[i];
      for(unsigned c(0);c<n_cats;++c){
        is >> _data.group_cat_post_prob[i][c];
      }
    }

    load_advance_gcp(is,end_label);

  } else if(dummy == section_id) {
    is >> dummy;
    if(dummy != end_label){
      load_die(is,end_label,dummy);
    }

  } else {
    die("unexpected input in sml model file");

  }

  register_param();

#ifdef DEBUG
  {
    const unsigned ps(param_size());
    simple_array<smlfloat> pcopy(ps);
    this->param_state(pcopy.begin());
    for(unsigned i(0);i<ps;++i){
      log_os << pcopy[i] << "\n";
      if( is_float_invalid(pcopy[i]) ){
        log_os << "ERROR:: subs_ml_model invalid input param: " << i << " " << pcopy[i] << "\n";
        abort();
      }
    }
  }
#endif
}



/// dtor is being explicitly declared here so that the auto_ptr's to forward
/// declared classes in the header can properly clean themselves up
///
subs_ml_model::
~subs_ml_model() {}



void
subs_ml_model::
reset_param(const PARAM_INIT_TYPE::index_t pinit){
  // everything else below depends on cat manager values... so this needs to be reset first:
  get_cat_manager_nonconst().reset_param(pinit);

  // rate(bg_gtor) and root depend on obs, so the order here is important!
  obs_nonconst().reset_param(pinit);

  get_time_gtor_nonconst().reset_param(pinit);
  get_bg_gtor_nonconst().reset_param(pinit);
  get_rate_gtor_nonconst().reset_param(pinit);
  get_root_gtor_nonconst().reset_param(pinit);
}



void
subs_ml_model::
report_time(std::ostream& os) const {

  const time_gtor& tg(get_time_gtor());
  const unsigned n_time_cats(tg.time_cat_size());
  const unsigned n_branches(tree().branch_size());
  std::vector<irv_t<smlfloat> > branch_time(n_branches);
  for(unsigned tc(0);tc<n_time_cats;++tc){
    for(unsigned i(0);i<n_branches;++i){
      branch_time[i] = tg.time_cat_branch_time_variance(i,tc);
    }

    os << "time_cat: " << tc << "\n";
    report_time_instance("indy_cat_bg_nuc_normalized_branch_time",branch_time,tree(),os);
    os << "\n";
  }
}



void
subs_ml_model::
report(std::ostream& os) const {

  const unsigned ps(param_size(PARAM_VIEW::INDY_MIN_PARAM));
  os << "Free parameters: " << ps << "\n";
  os << "State Size: " << state_size() << "\n";
  os << "\n\n";

  get_cat_manager().report(os);
  get_root_gtor().report(tree().root()->label().c_str(),os);
  report_time(os);
  get_rate_gtor().report(os);
}



void
subs_ml_model::
set_obs_cat_seq_state_counts(smlfloat const * const * const * c){
  obs_nonconst().set_obs_cat_seq_state_counts(c);
}



void
subs_ml_model::
set_obs_cat_seq_state_distro(prob_t const * const * const * c){
  obs_nonconst().set_obs_cat_seq_state_distro(c);
}



const prob_t * const *
subs_ml_model::
obs_seq_state_distro_cat(const unsigned cat_no) const {
  return obs_seq_state_distro_obs_cat(get_cat_manager().typed_cat_no_from_cat_no(cat_no,CAT_PARAM_TYPE::OBS));
}



static
void
check_obs_cat_no(const cat_manager& cm,
                 const unsigned obs_cat_no){
  const unsigned n_adsets(cm.assigned_data_set_size());
  const unsigned n_obs_cats(cm.typed_cat_size(CAT_PARAM_TYPE::OBS));
  if(obs_cat_no>n_obs_cats){
    die("invalid obs cat");
  }
  if(n_adsets>1 && n_adsets != n_obs_cats){
    die("Can't currently handle nontrivial obs cat+assigned data set combinations.");
  }
}



const prob_t * const *
subs_ml_model::
obs_seq_state_distro_obs_cat(const unsigned obs_cat_no) const {
  check_obs_cat_no(get_cat_manager(),obs_cat_no);
  return obs().seq_state_distro_obs_cat()[obs_cat_no];
}



const prob_t *
subs_ml_model::
obs_state_distro_cat(const unsigned cat_no) const {
  return obs_state_distro_obs_cat(get_cat_manager().typed_cat_no_from_cat_no(cat_no,CAT_PARAM_TYPE::OBS));
}



const prob_t *
subs_ml_model::
obs_state_distro_obs_cat(const unsigned obs_cat_no) const {
  check_obs_cat_no(get_cat_manager(),obs_cat_no);
  return obs().state_distro_obs_cat()[obs_cat_no];
}



void
subs_ml_model::
post_write_param_state()  {
  get_rate_gtor_nonconst().fix_input_param(); // moved up here b/c this relies on cat-manager update
  get_root_gtor_nonconst().parent_model_update();
}



smlfloat
subs_ml_model::
param_penalty() const {
  return get_rate_gtor().param_penalty()+get_time_gtor().param_penalty();
}



void
subs_ml_model::
register_param(){
  clear_register();
  base_t::register_param(get_time_gtor_nonconst());
  base_t::register_param(get_rate_gtor_nonconst());
  base_t::register_param(get_bg_gtor_nonconst());
  base_t::register_param(get_root_gtor_nonconst());
  base_t::register_param(get_cat_manager_nonconst());
}



void
subs_ml_model::
set_rate_gtor_model(const subs_ml_model_init_options& sml_opt){
  _rate_gtor.reset(rate_gtor_factory(rate_gtor_sml_share(*this),sml_opt.ropt,sml_opt.nopt,sml_opt.copt,sml_opt.c5t));
}



unsigned
subs_ml_model::
state_size() const { return RATE_GTOR_MODEL::state_size(rate_gtor_model()); }



unsigned
subs_ml_model::
submodel_size() const { return get_rate_gtor().submodel_size(); }



unsigned
subs_ml_model::
submodel_state_size(unsigned s) const { return get_rate_gtor().submodel_state_size(s); }



const notifier&
subs_ml_model::
obs_info_notifier() const { return obs(); }



const notifier&
subs_ml_model::
bg_gtor_notifier() const { return get_bg_gtor(); }

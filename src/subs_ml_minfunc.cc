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

// $Id: subs_ml_minfunc.cc 1179 2008-03-26 00:23:59Z ctsa $

/// \file

#include "cat_manager.h"
#include "condition_func.h"
#include "lhood_model.h"
//#include "lhood_model_em.h"
#include "lhood_model_prep.h"
#include "lhood_model_root.h"
#include "rate_gtor.h"
#include "root_gtor.h"
#include "site_data_fastlup.h"
#include "site_prob_maker.h"
#include "substk_exception.h"
#include "subs_ml_minfunc.h"
#include "subs_ml_model_min_options.h"
#include "subs_ml_print_util.h"
#include "time_gtor.h"
#include "util/general/log.h"
#include "util/math/array_util.h"
#include "util/math/test_float.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>


#ifdef DEBUG
#define MINFUNC_DEBUG
#endif


unsigned min_call_count(0);
unsigned dmin_call_count(0);
bool write_ping(true);

double last_lnp(std::log(0.)); ///< used during profiling



// occasional tasks during minimization
static
void
xminfunc_count_check(const smlfloat lnp,
                     const smlfloat,
                     const subs_ml_model& mdl,
                     const bool is_mstep = false){

  // return neg lnp for minimization
  static const unsigned iter_report_period(50);
  const unsigned counter(min_call_count+dmin_call_count);
  if ( counter%iter_report_period == 0 ){
    log_os << "\n";
    if(is_mstep) log_os << "mstep ";
    log_os << "lnp: "
           << std::setprecision(SUBS_ML_PRINT_LNP_PRECISION)
           << -lnp << "\n";
  }

  // backup model state
  const unsigned ping_pong_rate(180000/mdl.state_size());
  const unsigned ping_pong_first_jump(5*ping_pong_rate/6);
  if ( mdl.opt().is_pingpong &&
       (counter+ping_pong_first_jump)%ping_pong_rate == 0 ){
    std::string nfile(mdl.opt().pingpongfile);
    if( write_ping ) nfile += ".ping";
    else             nfile += ".pong";

    std::ofstream fos(nfile.c_str());
    mdl.store_state(fos);
    write_ping = ! write_ping;
  }
}



bool
subs_ml_minfunc_base::
is_val_computable(const smlfloat* p){

  min_p_to_param(p);
  return is_current_val_computable();
}



bool
subs_ml_minfunc_base::
is_current_val_computable() const {
  {
    const bool is_penalty(mdl().param_penalty() > 0.);
    if(is_penalty) return false;
  }

  {
    const unsigned n_cats(mdl().get_cat_manager().cat_size());

    bool is_ztime(false);
    const time_gtor& tg(mdl().get_time_gtor());
    for(unsigned i(0);i<n_cats;++i){
      is_ztime = is_ztime || tg.is_zero_time(i);
    }
    if(is_ztime) return false;
  }

  return true;
}



smlfloat
subs_ml_minfunc_base::
dval(const smlfloat* p,smlfloat* df){

  const unsigned ndim(dim());
  const bool xfv(is_inner_step());

  // used during profiling
  if( mdl().opt().max_steps ){
    if( (min_call_count+dmin_call_count) > mdl().opt().max_steps ) {
      for(unsigned i(0);i<ndim;++i) df[i] = 0.;
      return -last_lnp;
    }
  }

  // load minimizer parameters:
  min_p_to_param(p);

  // get current lnp
  smlfloat lnp,lnp_norm;
  get_lnprob(mdl(),lnp,lnp_norm);

  if( mdl().opt().max_steps ) last_lnp = lnp_norm;

  min_call_count++;
  xminfunc_count_check(lnp,lnp_norm,mdl(),xfv);

  // estimate deriv at low resolution
  for(unsigned i(0);i<ndim;++i){
    df[i] = get_deriv_estimate(p,lnp_norm,i);
    min_call_count++;
    xminfunc_count_check(lnp,lnp_norm,mdl(),xfv);
  }

  { // get norm of deriv:
    static const smlfloat norm_thresh(1.e-5);
    const smlfloat norm(std::sqrt(array_dot(df,df,ndim)));
    if(norm<norm_thresh){
      // estimate centered deriv
      for(unsigned i(0);i<ndim;++i){
        df[i] = (df[i] + get_deriv_estimate(p,lnp_norm,i,false))/2.;
        min_call_count++;
        xminfunc_count_check(lnp,lnp_norm,mdl(),xfv);
      }
    }
  }

  for(unsigned i(0);i<ndim;++i){ df[i] *= -1.; }

  dmin_call_count++;
  xminfunc_count_check(lnp,lnp_norm,mdl(),xfv);

#ifdef MINFUNC_DEBUG
  for(unsigned i(0);i<ndim;++i) { log_os << "i,p[i]: " << i << " " << p[i] << "\n"; }
  for(unsigned i(0);i<ndim;++i) { log_os << "i,df[i]: " << i << " " << df[i] << "\n"; }
  log_os << "return lnp/lnp_norm " << -lnp << " " << -lnp_norm << "\n";
#endif

  return -lnp_norm;
}



smlfloat
subs_ml_minfunc_base::
val(const smlfloat* p) {

  if( mdl().opt().max_steps ){
    if( (min_call_count+dmin_call_count) > mdl().opt().max_steps ) {
      return -last_lnp;
    }
  }

  min_p_to_param(p);

  // recalculate probability based on model parameters
  smlfloat lnp,lnp_norm;
  get_lnprob(mdl(),lnp,lnp_norm);

#ifdef MINFUNC_DEBUG
  {
    log_os << "minfunc:\n";
    const unsigned n(dim());
    for(unsigned i(0);i<n;++i) { log_os << "i,p[i]: " << i << " " << p[i] << "\n"; }
    log_os << "lnp/lnp_norm: " << -lnp << " " << -lnp_norm << "\n";
  }
#endif

  if( mdl().opt().max_steps ) last_lnp = lnp_norm;
  min_call_count++;
  xminfunc_count_check(lnp,lnp_norm,mdl(),is_inner_step());

  return -lnp_norm;
}



smlfloat
subs_ml_minfunc_base::
get_deriv_estimate(const smlfloat* p,
                   const smlfloat lnp1_norm,
                   const unsigned param_no,
                   const bool is_plus) const {

  static const smlfloat delta(1e-5);
  static const smlfloat deltamin(1e-7);

  const unsigned n(dim());
  simple_array<smlfloat> pcopy(n);
  std::copy(p,p+n,pcopy.ptr());

  const smlfloat delchunk(std::max(std::fabs(pcopy[param_no]*delta),deltamin));
  const smlfloat start_param(pcopy[param_no]);

  if(is_plus) pcopy[param_no] = start_param+delchunk;
  else        pcopy[param_no] = start_param-delchunk;

  //  log_os << "i,delchunk,pstart,pdel: " << param_no << " " << delchunk << " " << start_param << " " << pcopy[param_no] << "\n";

  subs_ml_model mdl_copy(mdl());
  min_p_to_param_x(pcopy.ptr(),po_from_mdl(mdl_copy));

  //  mdl_copy.param_state(pcopy,pv);

  //  log_os << "pdel_actual: " << pcopy[param_no] << "\n";

  smlfloat lnp2,lnp2_norm;
  get_lnprob(mdl_copy,lnp2,lnp2_norm);

  //  log_os << "lnp1_norm,lnp2_norm: " << lnp1_norm << " " << lnp2_norm << "\n";
  //  log_os << "returning: " << (lnp2_norm-lnp1_norm)/delchunk << "\n";

#if 0
  if(is_nonnegative_reflection){
    if(start_param<delchunk){
      if(is_plus) {
        const smlfloat standard((lnp2_norm-lnp1_norm)/delchunk);
        const smlfloat reflected((lnp1_norm-lnp2_norm)/(delchunk+start_param));
        return (standard+reflected)/2.;
      } else {
        pass_away("get_deriv_estimate: inconsistent arguments");
      }
    }
  }
#endif

  if(is_plus) return (lnp2_norm-lnp1_norm)/delchunk;
  else        return (lnp1_norm-lnp2_norm)/delchunk;
}



void
subs_ml_minfunc_base::
check_p(const smlfloat* p) const {
  const unsigned n(dim());
  for(unsigned i(0);i<n;++i){
    if(is_float_invalid(p[i])){
      std::ostringstream oss;
      oss << "minfunc.check_p(): invalid minimizer input i,p: " << i << " " << p[i] << "\n";
      for(unsigned j(0);j<n;++j){
        oss << "j,p: " << j << " " << p[j] << "\n";
      }
      throw substk_exception(oss.str().c_str());
    }
  }
}



subs_ml_minfunc::
subs_ml_minfunc(const site_data_fastlup& sdf,
                subs_ml_model& m,
                const PARAM_VIEW::index_t pv)
  : subs_ml_minfunc_base(m,pv),
    _sdf(sdf) {
  _nld.reset(new lhood_model_prep(m,sdf));
}



// see note in subs_ml_minfunc.h
subs_ml_minfunc::~subs_ml_minfunc() {}



void
subs_ml_minfunc::
get_lnprob(const subs_ml_model& mdl_copy,
           smlfloat& lnp,
           smlfloat& lnp_norm) const {
  get_lnprob_from_param_prep(mdl_copy,_sdf,*_nld,site_prob_maker_full(_sdf),lnp,lnp_norm);
}


#if 0
void
subs_ml_minfunc_em::
get_lnprob(const subs_ml_model& mdl_copy,
           smlfloat& lnp,
           smlfloat& lnp_norm) const {
  get_lnprob_from_param(mdl_copy,_edat,lnp,lnp_norm);
}
#endif



subs_ml_root_minfunc::~subs_ml_root_minfunc() {}



subs_ml_root_minfunc::
subs_ml_root_minfunc(const site_data_fastlup& sdf,
                     subs_ml_model& m,
                     const PARAM_VIEW::index_t pv)
  : subs_ml_minfunc_base(m,pv), _is_reset(false), _sdf(sdf),
    _root_spp(m.get_cat_manager().cat_size(),sdf.len,m.state_size()),
    _cf(m.get_cat_manager().cat_size()) {

  if(m.submodel_size()>1) die("subs_ml_root_minfunc(): invalid submodel min method");

 _nld.reset(new lhood_model_prep(m,sdf));

 const unsigned n_cats(m.get_cat_manager().cat_size());
 for(unsigned i(0);i<n_cats;++i){
   _cf[i].reset(m.get_rate_gtor().condition_func_factory());
   _cf[i]->data_init(m.tree(),sdf);
 }
}



void
subs_ml_root_minfunc::
reset(const subs_ml_model& m){
  get_root_node_site_partial_prob(m,_sdf,_root_spp,_cf.ptr());
  _is_reset=true;
}



void
subs_ml_root_minfunc::
get_lnprob(const subs_ml_model& mdl_copy,
           smlfloat& lnp,
           smlfloat& lnp_norm) const {

  if(! _is_reset) die("subs_ml_root_minfunc unformated...");

  const unsigned n_cats(mdl_copy.get_cat_manager().cat_size());
  for(unsigned c(0);c<n_cats;++c){
    _cf[c]->root_update(mdl_copy.get_root_gtor().cat_state_pdistro(c));
  }

  get_lnprob_from_param_prep(mdl_copy,_sdf,*_nld,
                             site_prob_maker_root(_sdf,_root_spp,_cf.ptr()),lnp,lnp_norm);
}



param_object_base&
subs_ml_root_minfunc::
po_from_mdl(subs_ml_model& m) const {
  return m.get_root_gtor_nonconst();
}



const param_object_base&
subs_ml_root_minfunc::
po_from_mdl(const subs_ml_model& m) const {
  return m.get_root_gtor();
}




subs_ml_rootcat_minfunc::~subs_ml_rootcat_minfunc() {}



subs_ml_rootcat_minfunc::
subs_ml_rootcat_minfunc(const site_data_fastlup& sdf,
                        subs_ml_model& m,
                        const PARAM_VIEW::index_t pv)
  : subs_ml_minfunc_base(m,pv), _is_reset(false), _cat_no(0),
    _sdf(sdf), _root_spp(sdf.len,m.state_size()),
    _othercat_site_prob(m.get_cat_manager().cat_size()-1,sdf.len),
    _cf(m.get_rate_gtor().condition_func_factory()) {

  if(m.submodel_size()>1) die("subs_ml_root_minfunc(): invalid submodel min method");

  _nld.reset(new lhood_model_prep(m,sdf));
  _cf->data_init(m.tree(),sdf);
}



void
subs_ml_rootcat_minfunc::
reset(const subs_ml_model& m,
      const unsigned cat_no){

  _cat_no=cat_no;

  // get partial probs @ root node for the category we're minimizing:
  //
  get_root_node_site_partial_prob_cat(m,_sdf,_cat_no,_root_spp.ptr(),*_cf);

  // get the full site probs for all the other categories:
  //
  const unsigned n_cats(m.get_cat_manager().cat_size());
  const unsigned n_sites(_sdf.len);
  site_prob_maker_full spm(_sdf);

  lhood_model_prep_cat_site_prob& csp(_nld->csp_prep);

  unsigned otherc(0);
  for(unsigned c(0);c<n_cats;++c){
    if(c==_cat_no) continue;
    spm.get_cat_site_prob(m,c,csp);
    for(unsigned i(0);i<n_sites;++i){
      _othercat_site_prob[otherc][i] = csp.site_prob[i];
    }
    otherc++;
  }

  _is_reset=true;
}



void
subs_ml_rootcat_minfunc::
get_lnprob(const subs_ml_model& mdl_copy,
           smlfloat& lnp,
           smlfloat& lnp_norm) const {

  if(! _is_reset) die("subs_ml_rootcat_minfunc unformated...");

  _cf->root_update(mdl_copy.get_root_gtor().cat_state_pdistro(_cat_no));

  site_prob_maker_rootcat spm(_sdf,_root_spp,_othercat_site_prob,_cat_no,*_cf);
  get_lnprob_from_param_prep(mdl_copy,_sdf,*_nld,spm,lnp,lnp_norm);
}



param_object_base&
subs_ml_rootcat_minfunc::
po_from_mdl(subs_ml_model& m) const {
  return m.get_root_gtor_nonconst().cat_param_object(_cat_no);
}



const param_object_base&
subs_ml_rootcat_minfunc::
po_from_mdl(const subs_ml_model& m) const {
  return m.get_root_gtor().cat_param_object(_cat_no);
}

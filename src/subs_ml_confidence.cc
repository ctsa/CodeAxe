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

// $Id: subs_ml_confidence.cc 1111 2008-01-26 01:11:52Z ctsa $

/// \file

#include "lhood_model.h"
#include "subs_ml_confidence.h"
#include "subs_ml_fisher.h"
#include "subs_ml_minimize.h"
#include "subs_ml_model_ci_info.h"
#include "subs_ml_model_min_options.h"
#include "substk_exception.h"
#include "util/general/die.h"
#include "util/general/lock_util.h"
#include "util/general/log.h"
#include "util/math/chi_sqr.h"
#include "util/math/matrix_invert.h"
#include "util/math/matrix_util_io.h"

#include <fstream>
#include <ostream>
#include <vector>


static
void
calc_param_ci(const site_data_fastlup& sdf,
              const subs_ml_model& mdl,
              const smlfloat lnp_conf_diff,
              const smlfloat start_lnp,
              const simple_array<bool>& init_is_train,
              const simple_array<smlfloat>& init_param,
              param_ci_info& pci){

  const unsigned ps(mdl.param_size());

  simple_array<bool> is_train(init_is_train);
  simple_array<smlfloat> param(init_param);

#define CI_DEBUG
#ifdef CI_DEBUG
  simple_array<smlfloat> param2(ps);

  log_os << "ps: " << ps << "\n";
#endif

  const unsigned i(pci.index);

  log_os << "calculating i: " << i << "\n";

  subs_ml_model mdl_ci(mdl);
  mdl_ci.opt().min = MIN::CONJ_DIR;

  // note that the lock change must come before changing param!, else
  // normalization could affect the changed param:
  is_train[i] = (! is_train[i]);
  mdl_ci.set_is_train_param_state(is_train.begin());

  /// /todo consider instead.. transform param trainable , +/-, then
  /// untransform to supply confidence interval in terms of real param
  /// value... or add public param min/max functions
  ///
  smlfloat param_delta;
  smlfloat new_param;

  // loop through increasing values of iter_delta to find a parameter
  // offset with a sufficiently large independent lnp diff before
  // searching for the dependent lnp diff
  //
  {
    static const smlfloat delta(1e-3);
    static const smlfloat min_param_delta(1e-5);
    static const smlfloat iter_delta_scale(1.5);

    const smlfloat ntol(mdl.opt().converge_value);

    for(smlfloat iter_delta(delta);iter_delta<100.;iter_delta *= iter_delta_scale) {
      param_delta = std::max(param[i]*iter_delta,min_param_delta*iter_delta);
      new_param = param[i] + param_delta;
      if((new_param<0 && param[i]>0) ||
         (new_param>0 && param[i]<0)){
        new_param = param[i] - param_delta;
      }

      std::swap(param[i],new_param);
      mdl_ci.set_param_state(param.begin());
      std::swap(param[i],new_param);

      smlfloat param_lnp,param_lnp_norm;
      get_lnprob_from_param(mdl_ci,sdf,param_lnp,param_lnp_norm);

      const smlfloat lnp_diff(param_lnp-start_lnp);

      if(-lnp_diff >= ntol*2000) break;
    }
  }

  try {
    log_os << "test param start " << new_param << "\n";

    mdl_ci.opt().is_extra_warnings=false; //make things quieter during iterative run
    sml_fullstep_min(sdf,mdl_ci,false);
    mdl_ci.opt().is_extra_warnings=true;

  } catch (substk_exception& s) {
    // if anything goes wrong, we want to know what the state of
    // mdl_ci is at the point of failure:
    //
    log_os << "FATAL:: SubsTK EXCEPTION:: " << s.what() << "\n";
    log_os << "...caught in " << __FILE__ << ":" << __LINE__ << "\n";
    log_os << "...dumping model:\n";
    mdl_ci.store_state(log_os);
    abort();
  }

  smlfloat param_lnp,param_lnp_norm;
  get_lnprob_from_param(mdl_ci,sdf,param_lnp,param_lnp_norm);

#ifdef CI_DEBUG
  mdl_ci.param_state(param2.begin());
  log_os << "test param finish " << param2[i] << "\n";
#endif

  pci.mean = param[i];
  pci.points.push_back(std::make_pair(new_param,param_lnp));


  // assuming normality, get variance of lnp wrt param...
  const smlfloat lnp_diff(param_lnp-start_lnp);


  if(lnp_diff>=0){
    // somthin's mussed up! ...add a more deliberate handler for this
    // situation later
    warning("invalid parameter response in ci calc");
  } else {

    const smlfloat variance(-(param_delta*param_delta)/(2.*lnp_diff));

    // then get significant param range
    const smlfloat target_delta_plus(std::sqrt(variance*(2.*lnp_conf_diff)));
    const smlfloat target_delta_minus(-target_delta_plus);

    pci.variance = variance;
    pci.ci_plus = target_delta_plus;
    pci.ci_minus = target_delta_minus;
  }

  log_os << "i,p,pdiff,lnp,lnpdiff,var: " << i
         << " " << param[i] << " " << new_param-param[i]
         << " " << start_lnp << " " << lnp_diff << " " << pci.variance
         << " +/-:  " << pci.ci_plus << " / "  << pci.ci_minus << "\n";

#ifdef CI_REALLY_DEBUG
  // nan test target_delta
  if(! (is_float_nan(pci.ci_plus) || is_float_nan(pci.ci_minus))){
    new_param=param[i]+pci.ci_plus;
    std::swap(param[i],new_param);
    mdl_ci.set_param_state(param.begin());
    std::swap(param[i],new_param);

    sml_fullstep_min(sdf,mdl_ci,false);
    smlfloat param_lnp_tplus,param_lnp_norm_tplus;
    get_lnprob_from_param(mdl_ci,sdf,param_lnp_tplus,param_lnp_norm_tplus);

    new_param=param[i]+pci.ci_minus;
    std::swap(param[i],new_param);
    mdl_ci.set_param_state(param.begin());
    std::swap(param[i],new_param);

    sml_fullstep_min(sdf,mdl_ci,false);
    smlfloat param_lnp_tminus,param_lnp_norm_tminus;
    get_lnprob_from_param(mdl_ci,sdf,param_lnp_tminus,param_lnp_norm_tminus);

    log_os << "start,plus,plusdiff,minus,minusdiff: "
           << start_lnp
           << " " << param_lnp_tplus  << " " << param_lnp_tplus-start_lnp
           << " " << param_lnp_tminus << " " << param_lnp_tminus-start_lnp << "\n";
  }
#endif
}



///
static
void
sml_conf_interval(const site_data_fastlup& sdf,
                  const subs_ml_model& mdl,
                  const smlfloat alpha,
                  const simple_array<bool>& ci_mask,
                  const bool is_single_param,
                  const std::string& outfile){

  assert(alpha < 1. && alpha > 0.);

  const smlfloat lnp_conf_diff(chi_sqr(1.-alpha,1.)/2.);

  // get lnp of original model, the model is assumed to be
  // the mle wrd data (sdf) for this calculation:
  //
  smlfloat start_lnp,start_lnp_norm;
  get_lnprob_from_param(mdl,sdf,start_lnp,start_lnp_norm);

  const unsigned ps(mdl.param_size());

  simple_array<bool> is_train(ps);
  mdl.is_train_param_state(is_train.begin());

  simple_array<smlfloat> param(ps);
  mdl.param_state(param.begin());

  // get pps:
  unsigned pps(0);
  for(unsigned i(0);i<ps;++i){
    if(! is_train[i]) continue;
    pps++;
  }

  for(unsigned i(0);i<ps;++i){

    if(! is_train[i]) continue;

    log_os << "checking i: " << i << "\n";

    param_ci_info pci;
    pci.index=i;
    pci.is_calculating=true;

    {
      file_locker fl(outfile.c_str());
      std::fstream fp(outfile.c_str(), std::fstream::in|std::fstream::out);

      // check that finished model file isn't already written:
      if(subs_ml_model_ci_info::check_state(fp)) break;

      fp.clear();
      fp.seekg(0,std::ios::beg);

      param_ci_info pci_tmp;

      bool is_bail(false);
      while(pci_tmp.load_state(fp)){
        if(pci_tmp.index == i){
          is_bail=true;
          break;
        }
      }
      if(is_bail) continue;

      fp.clear();
      pci.store_state(fp);
      fp << std::flush;
    }

    if(! ci_mask[i]){
      calc_param_ci(sdf,mdl,lnp_conf_diff,start_lnp,is_train,param,pci);
    }

    {
      file_locker fl(outfile.c_str());

      std::fstream fp(outfile.c_str(),std::fstream::in|std::fstream::out);
      fp.seekg(0, std::ios::end);
      pci.is_calculating = false;
      pci.store_state(fp);
      fp << std::flush;
    }

    if((! ci_mask[i]) && is_single_param) break;
  }


  {
    file_locker fl(outfile.c_str());
    std::fstream ifp(outfile.c_str(),std::fstream::in);

    // check that finished model file isn't already written:
    if(subs_ml_model_ci_info::check_state(ifp)) return;

    ifp.clear();
    ifp.seekg(0,std::ios::beg);

    unsigned linecount(0);
    while(ifp) { if(ifp.get() == '\n') linecount++; }
    if(linecount>=pps*2){
      subs_ml_model_ci_info cinfo(mdl.get_audit_info());
      cinfo.alpha = alpha;
      cinfo.base_lnp = start_lnp;
      ifp.clear();
      ifp.seekg(0,std::ios::beg);
      param_ci_info pci;
      while(pci.load_state(ifp)){
        if(! pci.is_calculating){
          cinfo.val.push_back(pci);
        }
      }
      std::sort(cinfo.val.begin(),cinfo.val.end());

      ifp.close();

      std::ofstream ofp(outfile.c_str(),std::fstream::out|std::fstream::trunc);
      cinfo.store_state(ofp);
    }
  }
}




void
get_subs_ml_model_conf_interval(const site_data_fastlup& sdf,
                                const subs_ml_model& mdl,
                                const conf_options& co,
                                const simple_array<bool>& ci_mask,
                                const std::string& outfile){

  // masking any free parameters from the fisher matrix is
  // meaningless... (according to ci_mask's intent)
  //
  if(co.is_fisher_method) {
#ifdef USE_LAPACK
    const unsigned ndim(mdl.param_size(PARAM_VIEW::TRAINABLE));
    std::vector<smlfloat> param_var(ndim);
    {
      smlfloat* fmat(new smlfloat[ndim*ndim]);

      get_fisher_matrix(sdf,mdl,fmat);

#ifdef DEBUG
      matrix_dump(fmat,ndim,log_os);
#endif

      matrix_invert(fmat,ndim);

#ifdef DEBUG
      matrix_dump(fmat,ndim,log_os);
#endif

      for(unsigned i(0);i<ndim;++i){ param_var[i] = fmat[i*(ndim+1)]; }
      delete [] fmat;
    }

    {
      assert(co.alpha < 1. && co.alpha > 0.);
      const smlfloat lnp_conf_diff(chi_sqr(1.-co.alpha,1.)/2.);

      subs_ml_model_ci_info cinfo;
      cinfo.alpha = co.alpha;

      const unsigned ps(mdl.param_size());

      simple_array<smlfloat> param(ps);
      mdl.param_state(param.begin());

      simple_array<bool> is_train(ps);
      mdl.is_train_param_state(is_train.begin());

      for(unsigned i(0),ti(0);i<ps;++i){
        if(! is_train[i]) continue;

        const smlfloat target_delta_plus(std::sqrt(param_var[ti]*(2.*lnp_conf_diff)));
        const smlfloat target_delta_minus(-target_delta_plus);

        param_ci_info pci;
        pci.index = i;
        pci.mean = param[i];
        pci.variance = param_var[ti];
        pci.ci_plus = target_delta_plus;
        pci.ci_minus = target_delta_minus;

        cinfo.val.push_back(pci);
        ti++;
      }

      std::ofstream fos(outfile.c_str());
      cinfo.store_state(fos);
    }
#else
    die("ci-fisher option requires lapack build");
#endif
  } else {
    sml_conf_interval(sdf,mdl,co.alpha,ci_mask,co.is_single_param,outfile);
  }
}

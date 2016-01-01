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

// $Id: subs_ml_minimize.cc 1185 2008-03-27 21:16:31Z ctsa $

/// \file

#include "cat_manager.h"
#include "lhood_model.h"
#include "lhood_model_em.h"
#include "root_gtor.h"
#include "subs_ml_minfunc.h"
#include "subs_ml_minimize.h"
#include "subs_ml_model_min_options.h"
#include "subs_ml_print_util.h"
#include "util/general/die.h"
#include "util/general/log.h"
#include "util/math/minimize_conj_direction.h"
#include "util/math/minimize_conj_gradient.h"
#ifdef USE_PRAXIS
#include "util/math/minimize_praxis.h"
#endif

#include <sstream>



static
void
run_conj_dir_minimizer(subs_ml_minfunc_base& mf,
                       const smlfloat start_ntol,
                       const smlfloat ntol,
                       const unsigned max_iter=0){

  static const smlfloat line_tol(1e-7);
  static const smlfloat start_ratio(0.05);
  static const smlfloat min_start_dist(1e-6);

  const unsigned ndim(mf.dim());
  const unsigned reset_iter(std::max(ndim/4,static_cast<unsigned>(10)));

  simple_array<smlfloat> p(ndim);
  simple_array<smlfloat> conj_dir(ndim*ndim);

  smlfloat iter_ntol(start_ntol);
  unsigned total_iter(0);

  while(true){
    mf.param_to_min_p(p.ptr());
    std::fill(conj_dir.begin(),conj_dir.end(),0.);
    for(unsigned i(0);i<ndim;++i) {
      const smlfloat start_dist( std::max(std::fabs(p[i]*start_ratio),min_start_dist) );
      conj_dir[i*(ndim+1)] = start_dist;
    }

    unsigned cd_iter(reset_iter);
    if(max_iter!=0) cd_iter=std::min(reset_iter,(max_iter-total_iter));

    smlfloat fmin;
    unsigned iter;
    smlfloat final_iter_delta_f(0.);

    minimize_conj_direction(p.ptr(),conj_dir.ptr(),mf,iter_ntol,ntol,line_tol,
                            fmin,iter,final_iter_delta_f,cd_iter);
    mf.min_p_to_param(p.ptr());

    total_iter += iter;
    log_os << "CONJ DIR MINIMIZER CYCLE/TOTAL ITERATIONS: " << iter << " " << total_iter << "\n";

    if(final_iter_delta_f <= ntol) break;
    if(max_iter != 0 && total_iter >= max_iter) {
      log_os << "CONJ DIR MINIMIZER REACHED MAX ITER LIMIT\n";
      break;
    }
  }
}



static
void
run_conj_grad_minimizer(subs_ml_minfunc_base& mf,
                        const smlfloat ntol,
                        unsigned reset_iter=0,
                        const unsigned max_iter=1000){

  if(reset_iter==0) reset_iter=max_iter;

  const unsigned ndim(mf.dim());
  simple_array<smlfloat> p(ndim);

  unsigned total_iter(0);

  while(true){
    if(total_iter+reset_iter>max_iter) reset_iter=max_iter-total_iter;

    mf.param_to_min_p(p.ptr());

    smlfloat fmin;
    unsigned iter;
    smlfloat final_iter_delta_f(0.);

    minimize_conj_gradient(p.ptr(),mf,ntol,fmin,iter,final_iter_delta_f,reset_iter);
    mf.min_p_to_param(p.ptr());

    total_iter += iter;

    log_os << "CONJ GRAD MINIMIZER CYCLE/TOTAL ITERATIONS: " << iter << " " << total_iter << "\n";

    if(final_iter_delta_f <= ntol) break;
    if(max_iter != 0 && total_iter >= max_iter) {
      log_os << "CONJ GRAD MINIMIZER REACHED MAX ITER LIMIT\n";
      break;
    }
  }
}



void
sml_fullstep_min_tol(const site_data_fastlup& sdf,
                     subs_ml_model& mdl,
                     const smlfloat start_ntol,
                     const smlfloat end_ntol,
                     const bool is_skip_hybrid_cg){

  const smlfloat norm(get_lnp_norm(mdl,sdf));
  const smlfloat end_ntol_norm(end_ntol/norm);
  const smlfloat start_ntol_norm(start_ntol/norm);

  min_call_count = 0;
  dmin_call_count = 0;

  subs_ml_minfunc mff(sdf,mdl,PARAM_VIEW::MIN_PARAM);
  subs_ml_minfunc mfg(sdf,mdl,PARAM_VIEW::INDY_MIN_PARAM);

  if (mdl.opt().min == MIN::CG_PRAXIS_HYBRID){
    // start out with a short conj gradient min
    static const unsigned max_iter(DEFAULT_CONJ_GRAD_START_ITER);
    run_conj_grad_minimizer(mfg,end_ntol_norm,max_iter,max_iter);

    const unsigned ndim(mff.dim());
    smlfloat* p(new smlfloat[ndim]);

    mff.param_to_min_p(p);
#ifdef USE_PRAXIS
    minimize_praxis(p,mff,end_ntol_norm,1.0,0.05);
#else
    pass_away("praxis minimization unavailable in this build");
#endif
    mff.min_p_to_param(p);

    delete [] p;

  } else if     (mdl.opt().min == MIN::CONJ_DIR){

    run_conj_dir_minimizer(mff,start_ntol_norm,end_ntol_norm);

  } else if (mdl.opt().min == MIN::CONJ_GRAD){

    static const unsigned reset_iter(60);
    run_conj_grad_minimizer(mfg,end_ntol_norm,reset_iter);

  } else if(mdl.opt().min == MIN::CG_CD_HYBRID){

    if(! is_skip_hybrid_cg){
      static const unsigned max_iter(DEFAULT_CONJ_GRAD_START_ITER);
      run_conj_grad_minimizer(mfg,end_ntol_norm,max_iter,max_iter);
    }

    run_conj_dir_minimizer(mff,start_ntol_norm,end_ntol_norm);
  } else {
    die("unknown minimization option in fullstep min");
  }

  log_os << "TOTAL MIN  CALLS: " << min_call_count <<"\n"
         << "TOTAL DMIN CALLS: " << dmin_call_count <<"\n";
}



void
sml_fullstep_min(const site_data_fastlup& sdf,
                 subs_ml_model& mdl,
                 const smlfloat start_lnp,
                 const bool is_start_tol){

  const smlfloat end_ntol(mdl.opt().converge_value);
  const smlfloat start_ntol(is_start_tol ? std::max((-start_lnp)*1e-3,end_ntol) : end_ntol);

  sml_fullstep_min_tol(sdf,mdl,start_ntol,end_ntol);
}



void
sml_root_cycle_min(const site_data_fastlup& sdf,
                   subs_ml_model& mdl,
                   const smlfloat start_lnp){

  const smlfloat end_ntol(mdl.opt().converge_value);
  const smlfloat start_ntol(std::max((-start_lnp)*1e-3,end_ntol));

  const smlfloat norm(get_lnp_norm(mdl,sdf));
  const smlfloat end_ntol_norm(end_ntol/norm);
  const smlfloat start_ntol_norm(start_ntol/norm);

  min_call_count = 0;
  dmin_call_count = 0;

#ifdef USE_ROOTCAT_MIN
  const unsigned n_cats(mdl.get_cat_manager().cat_size());
#endif

  // non-root step -> swap root param lock settings for full lock:
  // ...very ugly, but just expirmenting...
  //
  root_gtor& rg(mdl.get_root_gtor_nonconst());

  const unsigned ps(rg.param_size());
  simple_array<bool> is_train(ps);
  rg.is_train_param_state(is_train.begin());

  if (mdl.opt().min == MIN::CG_CD_HYBRID){

    smlfloat lnp,lnp_norm;

    {
      static const unsigned max_iter(DEFAULT_CONJ_GRAD_START_ITER);

      subs_ml_minfunc mfg(sdf,mdl,PARAM_VIEW::INDY_MIN_PARAM);
#ifdef USE_ROOTCAT_MIN
      subs_ml_rootcat_minfunc mfgr(sdf,mdl,PARAM_VIEW::INDY_MIN_PARAM);
#else
      subs_ml_root_minfunc mfgr(sdf,mdl,PARAM_VIEW::INDY_MIN_PARAM);
#endif

      for(unsigned i(0);i<2;++i){
        get_lnprob_from_param(mdl,sdf,lnp,lnp_norm);
        log_os << "PRE nroot-CG lnp: ";
        lnp_line_report(lnp,lnp_norm,sdf.total_count,log_os);

        rg.set_is_train_param_state(false);
        run_conj_grad_minimizer(mfg,end_ntol_norm,max_iter,max_iter/2);

        get_lnprob_from_param(mdl,sdf,lnp,lnp_norm);
        log_os << "PRE root-CG lnp: ";
        lnp_line_report(lnp,lnp_norm,sdf.total_count,log_os);

        rg.set_is_train_param_state(is_train.begin());

#ifdef USE_ROOTCAT_MIN
        for(unsigned c(0);c<n_cats;++c){
          mfgr.reset(mdl,c);
          run_conj_grad_minimizer(mfgr,end_ntol_norm,max_iter,max_iter/2);
        }
#else
        mfgr.reset(mdl);
        run_conj_grad_minimizer(mfgr,end_ntol_norm,max_iter,max_iter/2);
#endif
      }
    }

    {
      subs_ml_minfunc mff(sdf,mdl,PARAM_VIEW::MIN_PARAM);
#ifdef USE_ROOTCAT_MIN
      subs_ml_rootcat_minfunc mffr(sdf,mdl,PARAM_VIEW::MIN_PARAM);
#else
      subs_ml_root_minfunc mffr(sdf,mdl,PARAM_VIEW::MIN_PARAM);
#endif

      smlfloat ntol_norm_loop(start_ntol_norm);
      unsigned iter(0);
      while(true){

        smlfloat lnp1,lnp_norm1;
        get_lnprob_from_param(mdl,sdf,lnp1,lnp_norm1);
        log_os << "PRE nroot-CD cyc: " << iter << "  lnp: ";
        lnp_line_report(lnp1,lnp_norm1,sdf.total_count,log_os);

        rg.set_is_train_param_state(false);
        run_conj_dir_minimizer(mff,ntol_norm_loop,ntol_norm_loop);

        get_lnprob_from_param(mdl,sdf,lnp,lnp_norm);
        log_os << "PRE root-CD cyc: " << iter << "  lnp: ";
        lnp_line_report(lnp,lnp_norm,sdf.total_count,log_os);

        rg.set_is_train_param_state(is_train.begin());

#ifdef USE_ROOTCAT_MIN
        for(unsigned c(0);c<n_cats;++c){
          mffr.reset(mdl,c);
          run_conj_dir_minimizer(mffr,ntol_norm_loop,ntol_norm_loop);
        }
#else
        mffr.reset(mdl);
        run_conj_dir_minimizer(mffr,ntol_norm_loop,ntol_norm_loop);
#endif

        smlfloat lnp2,lnp_norm2;
        get_lnprob_from_param(mdl,sdf,lnp2,lnp_norm2);

        if((lnp_norm2-lnp_norm1) <= ntol_norm_loop) {
          if(lnp_norm2<lnp_norm1) {
            std::ostringstream oss;
            oss << "sml_rootcycle_min(): iteration " << iter << " lnL decreased from: "
                << lnp1 << " to " << lnp2 << "\n";
            warning(oss.str().c_str());
          }

          if(ntol_norm_loop*0.999 <= end_ntol_norm) break;
          static const smlfloat dfactor(.1);
          ntol_norm_loop *= dfactor;

          /// loop valuue 'snaps' to end value as soon as it gets close,
          /// which prevents an inefficient near-equal final iteration
          if((ntol_norm_loop*0.9)<= end_ntol_norm) ntol_norm_loop=end_ntol_norm;
        }
        ++iter;
      }
    }
  } else {
    die("unknown minimization option in rootcycle min");
  }

  log_os << "TOTAL MIN  CALLS: " << min_call_count <<"\n"
         << "TOTAL DMIN CALLS: " << dmin_call_count <<"\n";
}


#if 0
void
sml_mstep_min(const estep_dat_t& edat,
              subs_ml_model& mdl){

  const smlfloat norm(get_lnp_norm(mdl,edat));
  smlfloat ntol(mdl.opt().converge_value/norm);

  static const smlfloat mstep_tol_factor(1e2);

  ntol *= mstep_tol_factor; // looser tolerance b/c mstep is an "inner" iteration

  min_call_count = 0;
  dmin_call_count = 0;

  subs_ml_minfunc_em mfe(edat,mdl,PARAM_VIEW::MIN_PARAM);

  run_conj_dir_minimizer(mfe,ntol,ntol);

  log_os << "TOTAL MSTEP MIN  CALLS: " << min_call_count <<"\n"
         << "TOTAL MSTEP DMIN CALLS: " << dmin_call_count <<"\n";
}
#endif

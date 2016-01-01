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

// $Id: subs_ml_model_train.cc 1183 2008-03-27 02:12:11Z ctsa $

/// \file

#include "bi_tree.h"
#include "cat_manager.h"
#include "cat_post_prob.h"
#include "estep_dat_t.h"
#include "lhood_model.h"
#include "lhood_model_em.h"
#include "lhood_model_em_baz.h"
#include "lhood_model_em_pprob.h"
#include "subs_ml_model_train.h"
#include "subs_ml_minfunc.h"
#include "subs_ml_minimize.h"
#include "subs_ml_model_min_options.h"
#include "subs_ml_print_util.h"
#include "tree_probs.h"
#include "util/general/log.h"
#include "util/math/array_util.h"

#include <cstdlib>
#include <cmath>

#include <iomanip>
#include <ostream>


#if 0
///
static
void
sml_estep(const bi_tree_sml_type& tree,
          const tree_probs& tprob,
          const site_data_fastlup& sdf,
          estep_dat_t& edat,
          const unsigned N){

  update_etrans(tree,tprob,sdf,N,edat);
#if 0
   log_os << "eroot 1,2: " << edat.root[0] << " " << edat.root[1] << "\n";
   log_os << "etrans[0-4]:\n";
   for(unsigned i(0);i<4;++i){
     matrix_dump(edat.trans[i],N,log_os);
   }
#endif
}
#endif


// temporarily marking out em code:
#if 0

#ifdef DEBUG
static
void
check_edat(const estep_dat_t& edat,
           const smlfloat expect){

  smlfloat sum(0.);
  sum = array_sum(edat.root,edat.n_states);
  if(std::fabs(sum-expect)>0.001){
    log_os << "bad edat root sum:" << sum << "\n";
    abort();
  }

  for(unsigned i(0);i<edat.n_branches;++i){
    sum = array_sum(edat.trans[i],edat.n_states*edat.n_states);
    log_os << "i,sum,exp: " << i << " " << sum << " " << expect << "\n";
    if(std::fabs(sum-expect)>0.001){
      log_os << "bad edat branch sum: " << sum << " branch: " << i << "\n";
      abort();
    }
  }
}
#endif


static
void
sml_em_min(const site_data_fastlup& sdf,
           subs_ml_model& mdl,
           const smlfloat start_lnp){

  const unsigned n_states(mdl.state_size());
  const smlfloat ntol(mdl.opt().converge_value);

  const unsigned n_branches(mdl.tree().branch_size());

  estep_dat_t edat(n_states,n_branches);
  unsigned iter(0);

  smlfloat iter_lnp(start_lnp);

  while(true){
    iter++;
    {
      tree_probs tprob;
      tprob.init_probs(mdl);

      // control execution of temporary baz code:
#ifdef IS_BAZ_MODE
      update_etrans_baz(mdl,mdl.tree(),tprob,sdf,n_states);
      break;
#else
      update_etrans(mdl.tree(),tprob,sdf,n_states,edat);
#endif

      smlfloat mstep_lnp,mstep_lnp_norm;
      get_lnprob_from_param(mdl,edat,mstep_lnp,mstep_lnp_norm);
      log_os << "Mstep lnp: ";
      lnp_line_report(mstep_lnp,mstep_lnp_norm,sdf.total_count,log_os);

#ifdef DEBUG
      check_edat(edat,sdf.total_count);
#endif

#if 0
      mdl.rate_gtor().cf().adjust_edat(sdf);
#endif
      //lnp = get_lnprob_from_param(mdl,mdl.rate_gtor().cf().edat2);
      //log_os << "GREPGREP conditioned mstep lnp: " << std::setprecision(12) << lnp << "\n";
      //check_edat(mdl.rate_gtor().cf().edat2,n_states);
    }

    sml_mstep_min(edat,mdl);

    smlfloat new_lnp,new_lnp_norm;
    get_lnprob_from_param(mdl,sdf,new_lnp,new_lnp_norm);

    log_os << "EmIter" << iter << " lnp:  ";
    lnp_line_report(new_lnp,new_lnp_norm,sdf.total_count,log_os);

    if( mdl.opt().max_steps ){
      if( (min_call_count
           +dmin_call_count) > mdl.opt().max_steps ) {
        break;
      }
    }

    if (std::fabs(iter_lnp-new_lnp) <= ntol) break;

    iter_lnp = new_lnp;
  }
}
#endif


void
subs_ml_model::
clear_post_prob_data() {
  warning("Clearing model category posterior probs which are incompatible with data.");
  _data.group_cat_post_prob.clear();
  _data.group_label.clear();
}



void
subs_ml_model::
check_cat_post_prob_data_match(const site_data_fastlup& sdf){

  if(! is_gcp()) return;

  if(_data.group_cat_post_prob.dim1() != _data.group_label.size()){
    die("unexpected cat post prob info");
  }

  const unsigned n_groups_model(_data.group_label.size());
  const unsigned n_groups_data(sdf.n_groups);

  if(n_groups_model != n_groups_data){
    clear_post_prob_data();
    return;
  }

  for(unsigned g(0);g<n_groups_data;++g){
    if( _data.group_label[g] != sdf.group_label[g] ){
      clear_post_prob_data();
      return;
    }
  }
}



void
subs_ml_model::
update_model_obs_distro(const site_data_fastlup& sdf){

  const bool is_gcp_flag(is_gcp());

  const cat_manager& cm(get_cat_manager());
  const unsigned n_cats(cm.cat_size());
  const unsigned n_obs_cats(cm.typed_cat_size(CAT_PARAM_TYPE::OBS));

  simple_array<prob_t> cp_prior(n_cats);
  cm.cat_pdistro(cp_prior.ptr());

  const unsigned n_seqs(sdf.n_orgs);
  const unsigned n_states(state_size());
  const unsigned n_adsets(cm.assigned_data_set_size());
  simple_matrix3d<smlfloat> ssc(n_obs_cats,n_seqs,n_states,0);

  /// \todo ads_fix:
  ///
  if(n_adsets>1) die("can't handles multiple adsets");

  for(unsigned i(0);i<sdf.len;++i){
    const unsigned asize(sdf.data[i].group_assigned_data_set_count.size());
    for(unsigned j(0);j<asize;++j){
      const unsigned group_id(sdf.data[i].group_assigned_data_set_count[j].first.group_id);
      const unsigned adset_id(sdf.data[i].group_assigned_data_set_count[j].first.assigned_data_set_id);
      assert(adset_id==0);
      const unsigned count(sdf.data[i].group_assigned_data_set_count[j].second);

      const prob_t* cp(cp_prior.ptr());
      if(is_gcp_flag) cp=_data.group_cat_post_prob[group_id];

      for(unsigned c(0);c<n_cats;++c){
        const smlfloat pc(cp[c]*count);
        const unsigned obs_cat(cm.typed_cat_no_from_cat_no(c,CAT_PARAM_TYPE::OBS));
        for(unsigned t(0);t<n_seqs;++t) {
          assert( sdf.data[i].index[t] < n_states );
          // filter ambiguous states here:
          if(sdf.data[i].index[t] == n_states) continue;
          ssc[obs_cat][t][sdf.data[i].index[t]] += pc;
        }
      }
    }
  }

  set_obs_cat_seq_state_counts(ssc.ptr());
}



bool
subs_ml_model::
update_model_cat_expectation(const site_data_fastlup& sdf,
                             const smlfloat tol){

  const cat_manager& cm(get_cat_manager());
  if( cm.assigned_data_set_size() > 1 ){
    pass_away("get_model_cat_expectation(): can't handle multiple assigned data sets");
  }

  const unsigned n_cats(cm.cat_size());
  {
    const unsigned n_group_cats(cm.group_cat_size());
    if(n_cats != n_group_cats) die("post-prob cat update requires group cats only (for now)");
  }

  // get group posterior prob:
  //
  const unsigned n_groups(sdf.n_groups);
  simple_init_matrix<prob_t> group_cat_post_prob_rollback;

  const smlfloat norm_tol(tol/get_lnp_norm(*this,sdf));

  smlfloat lnp,lnp_norm;
  smlfloat last_lnp_norm(0.);
  const CPPM::index_t pm(CPPM::GROUP);
  simple_init_matrix<prob_t> ppi;

  // Allowing an iterative obs post-prob update might lock the
  // model into an early local minima.
  //
  // For MAXML obs_opt mode:
  //
  // In any model with multiple obs distro's, the likelihood could be
  // lowered by the post-prob update. The likelihood may also get
  // worse when we are working with approximate likelihoods, or it may
  // get slightly worse under theoretically correct EM conditions (no
  // multiple obs or approx lhood) because of floating point precision
  // limitations.  Therefore, the last set of posterior probabilities
  // is allowed to "rollback" if the likelihood is not improved.
  //
  // For EQ obs_opt mode:
  //
  // Are goal is to reach and maintain equilibrium w.r.t. obs
  // posterior prob updates, no matter what the effect on the model
  // likelihood.
  //
  for(unsigned iter(1);true;iter++){

    get_lnprob_from_param(*this,sdf,lnp,lnp_norm,pm,&ppi);
    make_cat_post_prob(n_cats,sdf,pm,ppi.ptr());

    log_os << "UpdateObsIter" << iter << "  lnp: ";
    lnp_line_report(lnp,lnp_norm,sdf.total_count,log_os);

    if(iter==1){
      _data.group_label.resize(n_groups);
      for(unsigned i(0);i<n_groups;++i){
        _data.group_label[i]=sdf.group_label[i];
      }
    } else {
      if       (opt().obs_update_mode == OBS_UPDATE_MODE::MAXML){
        if((lnp_norm-last_lnp_norm) < norm_tol){
          if(last_lnp_norm>lnp_norm){

            _data.group_cat_post_prob=group_cat_post_prob_rollback;
            update_model_obs_distro(sdf);

            get_lnprob_from_param(*this,sdf,lnp,lnp_norm,pm,&ppi);
            make_cat_post_prob(n_cats,sdf,pm,ppi.ptr());

            log_os << "UpdateObsRollback  lnp: ";
            lnp_line_report(lnp,lnp_norm,sdf.total_count,log_os);
          }
          return (iter>2);
        }

      } else if(opt().obs_update_mode == OBS_UPDATE_MODE::EQUIL){
        if(std::abs(lnp_norm-last_lnp_norm) < norm_tol){ return true; }


      } else die("unknown obs state optimization mode");
    }

    last_lnp_norm=lnp_norm;

    // store postp values from this iteration in case we have to
    // rollback:
    //
    group_cat_post_prob_rollback=_data.group_cat_post_prob;
    _data.group_cat_post_prob=ppi;

    update_model_obs_distro(sdf);
  }
}



static
void
sml_cat_postp_cycle_min(const site_data_fastlup& sdf,
                        subs_ml_model& mdl,
                        const smlfloat start_lnp,
                        const bool is_cat_em){

  static const unsigned max_cat_postp_iter(50);
  static const smlfloat ntol_reduction_factor(1e-2);

  const smlfloat ntol(mdl.opt().converge_value);
  const smlfloat iter1_ntol(-start_lnp*1e-4);

  smlfloat iter_ntol(iter1_ntol);

  smlfloat last_lnp(start_lnp);
  bool is_final_tolerance(false);

  unsigned iter(0);
  while(iter<max_cat_postp_iter){
    iter++;
    log_os << "\nCat Postp Iter: " << iter << " start\n\n";

    {
      bool is_updated(true);
      if(iter>1) {
        // these updates are relatively efficient, so switch the behavior to
        // always converge the update at the final tolerance. previous
        // mechanism used update_ntol commented out below:
        //        smlfloat update_ntol(iter_ntol);
        //        if(update_ntol*0.9 < ntol) update_ntol=ntol;
        //
        is_updated=mdl.update_model_cat_expectation(sdf,ntol);
      }

      if((! is_updated) && is_final_tolerance) break;
    }

    smlfloat end_ntol(iter_ntol);
    if(end_ntol*0.9 < ntol) {
      end_ntol=ntol;
      is_final_tolerance=true;
    }
    const smlfloat start_ntol(std::min(iter1_ntol,end_ntol*1e3));

    log_os << "Starting Cat Postp Iter Min: iter1/final/start/end tol: "
           << iter1_ntol << " " << ntol << " " << start_ntol << " " << end_ntol << "\n";

    if(is_cat_em) mdl.is_cat_mstep() = true;

    const bool is_skip_hybrid_cg(iter>1);

    sml_fullstep_min_tol(sdf,mdl,start_ntol,end_ntol,is_skip_hybrid_cg);

    if(is_cat_em) mdl.is_cat_mstep() = false;

    smlfloat lnp,lnp_norm;
    get_lnprob_from_param(mdl,sdf,lnp,lnp_norm);

    log_os << "Cat Postp Iter: " << iter << " Min finished lnp: ";
    lnp_line_report(lnp,lnp_norm,sdf.total_count,log_os);

    if ((lnp-last_lnp) <= ntol) break;

    last_lnp = lnp;
    if(! is_final_tolerance) iter_ntol *= ntol_reduction_factor;
  }
}



void
train_subs_ml_model(const site_data_fastlup& sdf,
                    subs_ml_model& mdl,
                    bool is_refine){

  // get initial lnp
  //
  smlfloat start_lnp,start_lnp_norm;
  get_lnprob_from_param(mdl,sdf,start_lnp,start_lnp_norm);

  const unsigned site_count(sdf.total_count);
  log_os << "Initial lnp: ";
  lnp_line_report(start_lnp,start_lnp_norm,site_count,log_os);

  // full max lhood minimization:
  //
  mdl.opt().is_extra_warnings=false; //make things quieter during iterative run

  const bool is_missing_obsp(mdl.get_cat_manager().is_missing_obs_partition());

  // mutually-exclusive non-trivial minimizations:
  //
  if       (mdl.opt().is_cat_em){
    if(mdl.opt().min == MIN::EM ||
       mdl.opt().is_rootcycle ||
       is_refine) die("incompatable minimization options");
    sml_cat_postp_cycle_min(sdf,mdl,start_lnp,mdl.opt().is_cat_em);

  } else if(mdl.opt().min == MIN::EM){
    if(mdl.opt().is_rootcycle ||
       is_missing_obsp ||
       is_refine) die("incompatable minimization options");
    die("EM minimization currently disabled");
    //    sml_em_min(sdf,mdl,start_lnp);

  } else if(mdl.opt().is_rootcycle){
    if(is_refine || is_missing_obsp) die("incompatible minimzation options");
    sml_root_cycle_min(sdf,mdl,start_lnp);

  } else if(is_missing_obsp){
    if(is_refine) die("incompatible minimization_options");
    sml_cat_postp_cycle_min(sdf,mdl,start_lnp,mdl.opt().is_cat_em);

  } else {
    sml_fullstep_min(sdf,mdl,start_lnp,!(is_refine));

  }

  mdl.opt().is_extra_warnings=true;
  smlfloat final_lnp,final_lnp_norm;
  get_lnprob_from_param(mdl,sdf,final_lnp,final_lnp_norm);

  log_os << "Final lnp: ";
  lnp_line_report(final_lnp,final_lnp_norm,site_count,log_os);
}

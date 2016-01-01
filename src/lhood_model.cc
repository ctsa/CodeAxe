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

// $Id: lhood_model.cc 1222 2008-05-22 23:10:06Z ctsa $

/// \file

#include "cat_manager.h"
#include "lhood_model.h"
#include "lhood_model_prep.h"
#include "lhood_model_prep_group.h"
#include "rate_gtor.h"
#include "site_data_fastlup.h"
#include "site_prob_maker.h"
#include "subs_ml_model.h"
#include "subs_ml_model_min_options.h"
#include "substk_exception.h"
#include "tree_probs.h"
#include "util/general/die.h"
#include "util/general/log.h"
#include "util/math/array_util.h"
#include "util/math/test_float.h"


namespace PVL {
extern "C" {
#include "pr.c"
}
}



smlfloat
get_lnp_norm(const subs_ml_model& mdl,
             const site_data_fastlup& sdf){

  const unsigned datasize(sdf.total_count);
  return static_cast<smlfloat>(datasize)*
    std::log(static_cast<smlfloat>(mdl.get_rate_gtor().state_size_conditioned()));
}



typedef site_data_fastlup_cell::group_assigned_data_set_vector::const_iterator gavi;



// all tprobs transposed!
void
get_lnprob_from_param_prep(const subs_ml_model& mdl,
                          const site_data_fastlup& full_sdf,
                          lhood_model_prep& nlp,
                          const site_prob_maker& spm,
                          smlfloat& lnp,
                          smlfloat& lnp_norm,
                          const CPPM::index_t cat_post_prob_mode,
                          simple_init_matrix<prob_t>* ppip){

  const unsigned full_count_len(full_sdf.len);

  const cat_manager& cm(mdl.get_cat_manager());

  const unsigned n_cats(cm.cat_size());
  const unsigned n_group_cats(cm.group_cat_size());
  const unsigned n_adsets(cm.assigned_data_set_size());
  const unsigned n_groups(full_sdf.n_groups);

  const bool is_use_groups(n_group_cats>1);

#if 0
  /// \todo make post-prob mode find the joint: P(site,group|m),
  /// instead of forcing the user to pick one of the marginals
  ///
  ///
  if(is_use_groups && cat_post_prob_mode == CPPM::SITE){
    pass_away("Can't use group cats w/ site cat post-prob assignments");
  }
  {
    const bool is_use_sites(n_cats>1 && n_cats!=n_group_cats);
    if(is_use_sites && cat_post_prob_mode == CPPM::GROUP){
      pass_away("Can't use site cats w/ group cat post-prob assignments");
    }
  }

  // no efficiency constraint in cat post prob mode, so do what's safest/simplest/easiest:
  //
  simple_init_matrix<prob_t>& ppi(*ppip);
  if       (cat_post_prob_mode == CPPM::SITE){
    ppi.init(full_count_len,n_cats,0.);
  } else if(cat_post_prob_mode == CPPM::GROUP){
    ppi.init(n_groups,n_cats,0.);
  }
#endif

  // get category distros
  cm.cat_pdistro(nlp.cat_pdistro.begin());

  for(unsigned a(0);a<n_adsets;++a){
    cm.adset_cat_pdistro(nlp.adset_cat_pdistro[a],a);
  }

  const bool is_gppi(mdl.is_gcp());
  const prob_t* const * gppi(0);

  cm.group_cat_pdistro(nlp.group_cat_pdistro.begin());
  for(unsigned a(0);a<n_adsets;++a){
    cm.adset_group_cat_pdistro(nlp.adset_group_cat_pdistro[a],a);
  }

  { // marginalize over all adsets to get group_cat priors conditioned
    // on group-id:
    const prob_t* const * pag(nlp.prob_adset_on_group());
    for(unsigned g(0);g<n_groups;++g){
      for(unsigned gc(0);gc<n_group_cats;++gc){
        prob_t& gp(nlp.group_id_group_cat_pdistro[g][gc]);
        gp = 0.;
        for(unsigned a(0);a<n_adsets;++a){
          gp += nlp.adset_group_cat_pdistro[a][gc]*pag[g][a];
        }
        gp = std::log(gp);
      }
    }
  }

  group_data& gd(nlp.gd());
  std::fill(gd.group_prob.begin(),gd.group_prob.end(),PVL::pr_zero());
  if(mdl.opt().is_cat_em && is_gppi){ gppi=mdl.get_gcp(); }


  lnp = 0.;

  for(unsigned c(0);c<n_cats;++c){

    spm.get_cat_site_prob(mdl,c,nlp.csp_prep);

    const unsigned group_cat(cm.group_cat_no(c));

    const bool is_new_group(c==0 || group_cat != cm.group_cat_no(c-1));
    if(is_new_group){
      std::fill(gd.group_prob_tmp.begin(),gd.group_prob_tmp.end(),0.);

      for(unsigned a(0);a<n_adsets;++a){
        nlp.ads_group_cat_prior_norm[a] = 0.;
        for(unsigned i(0);i<full_count_len;++i){
          nlp.site_cat_mix_prob[a][i] = 0.;
        }
      }
    }

    const prob_t* site_prob(nlp.csp_prep.site_prob);

#if 0
    if(cat_post_prob_mode == CPPM::SITE){
      for(unsigned i(0);i<full_count_len;++i){
        ppi[i][c] += site_prob[i]*nlp.cat_pdistro[c];
      }
    }
#endif

    for(unsigned a(0);a<n_adsets;++a){
      if(nlp.csp_prep.is_adset_using_cat[a][c]){
        const prob_t& cat_prior(nlp.adset_cat_pdistro[a][c]);
        nlp.ads_group_cat_prior_norm[a] += cat_prior;
        for(unsigned i(0);i<full_count_len;++i){
          nlp.site_cat_mix_prob[a][i] += site_prob[i]*cat_prior;
        }
      }
    }


    const bool is_last_of_group_cat((c+1)==n_cats || group_cat != cm.group_cat_no(c+1));

    if(is_last_of_group_cat){
      // logify site-mix probs and normalize prior probs to sum to 1 over the group
      for(unsigned a(0);a<n_adsets;++a){
        array_log(nlp.site_cat_mix_prob[a],full_sdf.len);

        const smlfloat lagnorm(std::log(nlp.ads_group_cat_prior_norm[a]));
        for(unsigned i(0);i<full_sdf.len;++i){
          nlp.site_cat_mix_prob[a][i] -= lagnorm;
        }
      }

      for(unsigned i(0);i<full_sdf.len;++i){
        gavi g(full_sdf.data[i].group_assigned_data_set_count.begin());
        const gavi g_end(full_sdf.data[i].group_assigned_data_set_count.end());
        for(;g!=g_end;++g){
          const site_data_fastlup_code::group_id_type group_id(g->first.group_id);
          const site_data_fastlup_code::assigned_data_set_id_type adset_id(g->first.assigned_data_set_id);
          const unsigned count(g->second);
          gd.group_prob_tmp[group_id] += nlp.site_cat_mix_prob[adset_id][i]*static_cast<smlfloat>(count);
        }
      }

      if(mdl.is_cat_mstep() && is_gppi){
        for(unsigned g(0);g<n_groups;++g){
          const smlfloat lngcp(nlp.group_id_group_cat_pdistro[g][group_cat]);
          lnp += gppi[g][group_cat]*(gd.group_prob_tmp[g]+lngcp);
        }
      } else {
        if(is_use_groups){
          for(unsigned g(0);g<n_groups;++g){
            const smlfloat lngcp(nlp.group_id_group_cat_pdistro[g][group_cat]);
            gd.group_prob[g] = PVL::pr_add(gd.group_prob[g],PVL::pr_nats2pr(-(gd.group_prob_tmp[g]+lngcp)));
          }
        }
      }
#if 0
      if(cat_post_prob_mode == CPPM::GROUP){
        for(unsigned g(0);g<n_groups;++g){
          ppi[g][group_cat] = gd.group_prob_tmp[g]+lngcp;
        }
      }
#endif
    }
  }

  if(! (mdl.is_cat_mstep() && is_gppi) ){
    if(is_use_groups){
#if 0
      // this might be faster...
      PVL::pr_t prob(PVL::pr_unity());
      for(unsigned g(0);g<n_groups;++g){
        prob = PVL::pr_multiply(prob,gd.group_prob[g]);
      }
      lnp = -PVL::pr_pr2nats(prob);
#else
      for(unsigned g(0);g<n_groups;++g){
        lnp += -PVL::pr_pr2nats(gd.group_prob[g]);
      }
#endif
    } else {
      for(unsigned g(0);g<n_groups;++g){
        lnp += gd.group_prob_tmp[g];
      }
    }
  }

  const smlfloat norm(get_lnp_norm(mdl,full_sdf));

  // apply penalties for out-of-range parameters
  const smlfloat pen(mdl.param_penalty());
  if(pen>0.){
    log_os << "WARNING:: Applying parameter out-of-range penalty: " << pen << "\n";
    lnp -= pen*norm;
  }

  lnp_norm = lnp/norm;

  if( is_float_nan(lnp_norm)){
    log_os << "FATAL:: Invalid lnp_norm: " << lnp_norm << "\n";
    log_os << "Dumping model:\n";
    mdl.store_state(log_os);
    throw substk_exception("get_lnprob_from_param_prep(): Invalid lnp_norm");
  } else if( is_float_inf(lnp_norm)){
    log_os << "WARNING:: get_lnprob returning lnp_norm: " << lnp_norm << "\n";
    log_os << "          attempting to proceed...\n";
    mdl.store_state(log_os);
  } else if( lnp > 0. ){
    log_os << "VERY FATAL:: Positive lnp!!!: " << lnp << "\n";
    log_os << "Dumping model:\n";
    mdl.store_state(log_os);
    abort();
  }
}




void
get_lnprob_from_param(const subs_ml_model& mdl,
                     const site_data_fastlup& sdf,
                     smlfloat& lnp,
                     smlfloat& lnp_norm,
                     const CPPM::index_t cat_post_prob_mode,
                     simple_init_matrix<prob_t>* ppip){

  // parameter independent setup:
  lhood_model_prep nlp(mdl,sdf);

  site_prob_maker_full spm(sdf);

  get_lnprob_from_param_prep(mdl,sdf,nlp,spm,lnp,lnp_norm,cat_post_prob_mode,ppip);
}


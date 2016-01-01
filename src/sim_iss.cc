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

// $Id: sim_iss.cc 1193 2008-03-29 03:14:54Z ctsa $

/// \file

#include "cat_info.h"
#include "cat_manager.h"
#include "lhood_site.h"
#include "rate_gtor.h"
#include "sim_iss.h"
#include "sim_util.h"
#include "simple_util.h"
#ifdef DIRECT_SITE_CONV
#include "site_code.h"
#endif
#include "subs_ml_ptol.h"
#include "tree_probs.h"
#include "util/general/die.h"
#include "util/general/log.h"
#include "util/math/random_util.h"

#include <ostream>
#include <iomanip>
#include <vector>


/// \todo switch to an mcmc sampling scheme so that this method can be
/// applied to larger simulations
///
void
simulate_data_indy_site_sampler(const sim_options& sim_opt,
                                const subs_ml_model& mdl,
                                nuc_seq_data& nsd){

  static const unsigned MAX_SITES(1000000);

  const unsigned n_states(mdl.state_size());
  const cat_manager& cm(mdl.get_cat_manager());
  const unsigned n_cats(cm.cat_size());

  {
    const unsigned n_group_cats(cm.group_cat_size());
    if(n_group_cats>1) {
      die("Group categories not supported in iss simulator");
    }
  }

  std::vector<tree_probs> tprob(n_cats);

  for(unsigned c(0);c<n_cats;++c){
    tprob[c].init_probs(mdl,c);
  }

  const unsigned n_leaves(mdl.tree().leaf_size());

  unsigned n_sites(1);
  for(unsigned t(0);t<n_leaves;++t){
    /// \todo add overflow safeguards here
    n_sites *= n_states;
  }

  simple_array<prob_t> site_prob_cdf(n_sites*n_cats);

  const unsigned ping_size(n_sites/100);
  unsigned ping_perc(0);
  bool is_sim_ping(false);
  if(n_sites >= MAX_SITES) {
    is_sim_ping=true;
    log_os << "Starting sim_distro calculation:\n";
  }

  simple_array<prob_t> cp(n_cats);
  cm.cat_pdistro(cp.ptr());

  prob_t sum(0.);
  simple_array<site_code::index_type> leaf_id(n_leaves);
  for(unsigned c(0);c<n_cats;++c){
    const prob_t cat_prior(cp[c]);
    for(unsigned i(0);i<n_sites;++i){
      for(unsigned t(0),ii(i);t<n_leaves;++t){
        leaf_id[t] = ii%n_states;
        ii /= n_states;
      }

      const prob_t site_prob = get_site_prob_single_cat(leaf_id.ptr(),
                                                        mdl.tree(),
                                                        tprob[c],
                                                        n_states);

      sum += site_prob*cat_prior;
      site_prob_cdf[c*n_sites+i] = sum;

      if(is_sim_ping){
        if(((i+1)%ping_size)==0) {
          log_os << "cat: " << c << " sim_distro: " << (++ping_perc) << "% complete\n";
        }
      }
    }
  }


  log_os << "SITE PROB SUM= " << sum << "\n";
#ifdef DEBUG
  if(std::fabs(sum-1.)>SUBS_ML_PTOL){
    log_os << "Simulator prob distro sum deviates from unity: "
           << std::setprecision(12) << sum << "\n";
    abort();
  }
#endif

#ifdef DIRECT_SITE_CONV
  // write simulation out to sd:
  //
  sim_site_data_init(mdl.tree(),mdl.rate_gtor_model(),sd);

  static const char group_label[] = "sim_iss_group_0";
  sd.group_label.push_back(group_label);

  // sample from site_prob_cdf to produce simulated data, store
  // sim data in sd:
  for(unsigned i(0);i<sim_opt.size;++i){
    static const unsigned group_no(0);
    const unsigned ranindex(random_cdf_variate(site_prob_cdf.ptr(),n_sites*n_cats));
    const unsigned site_id(ranindex%n_sites);
    const unsigned cat_no(ranindex/n_sites);

    site_code sc(n_leaves);
    for(unsigned t(0),ii(site_id);t<n_leaves;++t){
      sc.set_taxid(t,ii%n_states);
      ii /= n_states;
    }

    unsigned observed_cat_no(UNASSIGNED_CAT_ID);
    if(sim_opt.assigned_cat_prob>0.){
      const prob_t r(random_uniform());
      if(r<sim_opt.assigned_cat_prob) observed_cat_no=cat_no+1;
    }

    sd.site_count[sc][count_data_code(group_no,observed_cat_no)]++;
  }
#else
  static const unsigned n_groups(1);

  sim_seq_data_init(mdl.tree(),n_groups,n_cats,"sim_iss",cm,nsd);

  const RATE_GTOR_MODEL::index_t rgm(mdl.rate_gtor_model());
  const SITE_MODEL::index_t smodel(RATE_GTOR_MODEL::convert_to_site_model(rgm));
  const unsigned bsize(SITE_MODEL::base_size(smodel));
  const unsigned nuc_size(sim_opt.size*bsize);

  align_dat_t& na(nsd.dat[0]);
  na.seq.init(n_leaves,nuc_size);
  na.info.resize(nuc_size);

  unsigned codon_pos[SITE_MODEL::MAX_BASE_SIZE];
  SITE_MODEL::codon_position(smodel,codon_pos);

  for(unsigned i(0);i<sim_opt.size;++i){
    const unsigned ranindex(random_cdf_variate(site_prob_cdf.ptr(),n_sites*n_cats));
    const unsigned site_id(ranindex%n_sites);
    const unsigned cat_no(ranindex/n_sites);

    // relies on assumption that sites are non-overlapping:
    for(unsigned t(0),ii(site_id);t<n_leaves;++t){
      SITE_MODEL::decode_nuc(smodel,ii%n_states,na.seq[t]+i*bsize);
      ii /= n_states;
    }

    unsigned observed_cat_no(UNASSIGNED_CAT_ID);
    if(sim_opt.assigned_cat_prob>0.){
      const prob_t r(random_uniform());
      if(r<sim_opt.assigned_cat_prob){ observed_cat_no=cat_no+1; }
    }

    for(unsigned j(0);j<bsize;++j){
      const unsigned bi((i*bsize)+j);

      align_col_t& nai(na.info[bi]);

      nai.class_no=observed_cat_no;
      nai.codon_pos=codon_pos[j];
      nai.is_continuous = (bi!=0);
    }
  }

#endif
}

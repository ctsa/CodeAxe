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

// $Id: sim_discrete.cc 1193 2008-03-29 03:14:54Z ctsa $

/// \file

#include "bi_tree.h"
#include "cat_info.h"
#include "cat_manager.h"
#include "rate_gtor.h"
#include "sim_discrete.h"
#include "sim_util.h"
#include "simple_util.h"
#ifdef DIRECT_SITE_CONV
#include "site_code.h"
#endif
#include "time_gtor.h"
#include "tree_probs.h"
#include "util/general/log.h"
#include "util/math/prob_util.h"
#include "util/math/random_util.h"

#include <ostream>
#include <sstream>
#include <vector>

using namespace std;



static
void
rate_to_discrete_time_anysub_prob(const smlfloat* rates,
                                  prob_t* anysub_prob,
                                  const smlfloat time_unit,
                                  const unsigned n_states){

  for(unsigned i(0);i<n_states;++i) anysub_prob[i] = -time_unit*rates[i*(1+n_states)];
}



static
void
rate_to_discrete_time_anysub_cdf(const smlfloat* rates,
                                 prob_t* anysub_cdf,
                                 const unsigned n_states){

  for(unsigned i(0);i<n_states;++i){
    for(unsigned j(0);j<n_states;++j){
      if(j==i) continue;
      anysub_cdf[j+(j<i?0:-1)+i*(n_states-1)] = rates[j+i*n_states];
    }
    pdistro_norm(anysub_cdf+i*(n_states-1),anysub_cdf+(i+1)*(n_states-1));
    pdistro_to_cdf_inplace(anysub_cdf+i*(n_states-1),n_states-1);
  }
}



static
void
setup_anysub_rates(prob_t** anysub_prob,
                   prob_t** anysub_cdf,
                   const smlfloat time_unit,
                   const subs_ml_model& mdl,
                   const unsigned branch_no){

  const unsigned n_states(mdl.state_size());
  const rate_gtor& rg(mdl.get_rate_gtor());
  const cat_manager& cm(mdl.get_cat_manager());
  const unsigned n_cats(cm.cat_size());

  const time_gtor& tg(mdl.get_time_gtor());
  const smlfloat etime(tg.cat_expected_branch_time(branch_no));

  smlfloat* rates(new smlfloat[n_states*n_states]);

  for(unsigned c(0);c<n_cats;++c){
    const unsigned bcs(cm.get_branch_cat_set(c,branch_no));

    rg.rates(rates,rates_func_options(rates_func_options_base(c),bcs));

    rate_to_discrete_time_anysub_prob(rates,anysub_prob[c],time_unit,n_states);
    rate_to_discrete_time_anysub_cdf(rates,anysub_cdf[c],n_states);

    const smlfloat cat_branch_scale(tg.branch_time(branch_no,c)/etime);
    for(unsigned i(0);i<n_states;++i) anysub_prob[c][i] *= cat_branch_scale;
  }
  delete [] rates;
}



static
void
rates_diag_mutonly(const rate_gtor& r,
                   const unsigned cat,
                   const unsigned branch_cat_set,
                   smlfloat* rates_diag){

  const unsigned n_states(r.state_size());

  smlfloat* rates(new smlfloat[n_states*n_states]);
  r.rates(rates,rates_func_options(rates_func_options_base(cat,true,true),branch_cat_set));

  for(unsigned i(0);i<n_states;++i) rates_diag[i] = rates[i*(1+n_states)];
  delete [] rates;
}



static
void
setup_anysub_neutral_rates(prob_t** anysub_prob_neutral,
                           const smlfloat time_unit,
                           const subs_ml_model& mdl,
                           const unsigned branch_no){

  const unsigned n_states(mdl.state_size());
  const rate_gtor& rg(mdl.get_rate_gtor());
  const cat_manager& cm(mdl.get_cat_manager());
  const unsigned n_cats(cm.cat_size());

  const time_gtor& tg(mdl.get_time_gtor());
  const smlfloat etime(tg.cat_expected_branch_time(branch_no));

  smlfloat* rates_diag(new smlfloat[n_states]);

  for(unsigned c(0);c<n_cats;++c){
    const unsigned branch_cat_set(cm.get_branch_cat_set(c,branch_no));

    rates_diag_mutonly(rg,c,branch_cat_set,rates_diag);

    for(unsigned i(0);i<n_states;++i) anysub_prob_neutral[c][i] = -time_unit*rates_diag[i];

    const smlfloat time_cat_branch_scale(tg.branch_time(branch_no,c)/etime);
    for(unsigned i(0);i<n_states;++i) anysub_prob_neutral[c][i] *= time_cat_branch_scale;
  }
  delete [] rates_diag;
}


struct sim_info_t {
  vector<unsigned> cat_seq;
  prob_t** anysub_prob;
  prob_t** anysub_cdf;
  smlfloat time_unit_normed;              ///< discrete time unit
  unsigned n_states;
};



static
void
simulate_steps(vector<unsigned>& endseq,
               const sim_info_t& si,
               const bool is_report_sim_time,
               const prob_t* const * anysub_prob_neutral,
               const unsigned nsteps,
               const smlfloat last_fracstep,
               std::vector<smlfloat>& cat_sim_time){

  const unsigned n_cats(cat_sim_time.size());
  for(unsigned i(0);i<n_cats;++i) cat_sim_time[i] = 0.;

  const unsigned lsize(endseq.size());
  smlfloat fracstep(1.);

  for(unsigned step(0);step<(nsteps+1);++step){
    if(step==nsteps) fracstep=last_fracstep;

    for(unsigned i(0);i<lsize;++i){
      // lookup the site and group rate categories:
      const unsigned c(si.cat_seq[i]);
      const double r(random_uniform());
      const unsigned tmpl(endseq[i]);

      // record number of attempted (neutral) mutations
      if(is_report_sim_time){
        if(r<(anysub_prob_neutral[c][tmpl]*fracstep)) cat_sim_time[c] += 1.;
      }

      if(r<(si.anysub_prob[c][tmpl]*fracstep)){
        unsigned tmpi = random_cdf_variate(si.anysub_cdf[c]+tmpl*(si.n_states-1),si.n_states-1);
        if(tmpl<=tmpi) tmpi += 1;
        endseq[i] = tmpi;
      }
    }
  }
}



/// \param anysub_prob_neutral prob of any substitution during the
/// discrete time unit, absent selection
///
static
void
simulate_discrete_time_branch(const vector<unsigned>& startseq,
                              vector<unsigned>& endseq,
                              const smlfloat time,
                              const sim_info_t& si,
                              const bool is_report_sim_time,
                              const prob_t* const * anysub_prob_neutral,
                              std::vector<smlfloat>& cat_sim_time){

  const smlfloat time_tmp = time/si.time_unit_normed;
  const unsigned nsteps = static_cast<unsigned>(time_tmp);
  const smlfloat fracstep = time_tmp-static_cast<smlfloat>(nsteps);

  log_os << "calc branch: steps,frac: " << nsteps << " " << fracstep << "\n";

  endseq = startseq;

  simulate_steps(endseq,si,is_report_sim_time,anysub_prob_neutral,nsteps,fracstep,cat_sim_time);
}



// non-context, discrete time sim
void
simulate_data_discrete(const sim_options& sim_opt,
                       const subs_ml_model& mdl,
                       nuc_seq_data& nsd){

  const unsigned n_states(mdl.state_size());
  const rate_gtor& rg(mdl.get_rate_gtor());
  const cat_manager& cm(mdl.get_cat_manager());
  const unsigned n_cats(cm.cat_size());

  sim_info_t sim_info;

  sim_info.n_states = n_states;

  const bi_tree& tree(mdl.tree());
  const unsigned n_nodes(tree.node_size());
  const unsigned n_branches(tree.branch_size());
  const unsigned n_leaves(tree.leaf_size());

   // assign sites to groups:
  vector<unsigned> group_seq;
  get_sim_group_seq(group_seq,sim_opt);

  const unsigned n_groups(group_seq.back()+1);

  // assign cats to sites
  get_sim_cat_seq(sim_info.cat_seq,group_seq,sim_opt.size,mdl);

  // get root sequence:
  vector<vector<unsigned> > treenode_seq(n_nodes);
  get_ancestral_seq(mdl,sim_opt.size,sim_info.cat_seq,treenode_seq[tree.root()->node_id()]);


  // convert instantaneous rate to the prob of exchange in the discrete time unit:
  //
  //  anysub_prob - at any one site, the probability of changing state
  //  in 1 time unit
  //
  //  anysub_cdf - at any one site, the cdf for the probability of a
  //  different state, given that a change occurs
  //

  /// \todo see if this should be divided by 3 for codons - although this leads to
  ///       a pretty minor difference, it would be more consistent, no?
  sim_info.time_unit_normed = SIMULATE_DISCRETE_TIME_UNIT;

  // setup neutral substitution count data structures:
  std::vector<std::vector<smlfloat> > cat_sim_time(n_branches);
  simple_init_matrix3d<prob_t> anysub_prob_neutral;
  if(sim_opt.is_report_time){
    anysub_prob_neutral.init(n_branches,n_cats,n_states);
    for(unsigned i(0);i<n_branches;++i){
      setup_anysub_neutral_rates(anysub_prob_neutral[i],sim_info.time_unit_normed,mdl,i);
    }
    for(unsigned i(0);i<n_branches;++i){
      cat_sim_time[i].resize(n_cats);
    }
  }

  {
    // anysub_prob; ///< prob of any substitution during the discrete time unit
    // anysub_cdf;  ///< cdf of substitution in discrete time unit
    simple_matrix3d<prob_t> anysub_prob(n_branches,n_cats,n_states);
    simple_matrix3d<prob_t> anysub_cdf(n_branches,n_cats,n_states*(n_states-1));

    for(unsigned i(0);i<n_branches;++i){
      setup_anysub_rates(anysub_prob[i],
                         anysub_cdf[i],
                         sim_info.time_unit_normed,mdl,i);
    }

    // 3. propagate down branches to produce present day sequences:
    for(unsigned b(0);b<n_branches;++b){
      const bi_tree_node* branch_node(tree.branch_node(b));
      const smlfloat etime(mdl.get_time_gtor().cat_expected_branch_time(b));
      const prob_t* const * apn_neutral(anysub_prob_neutral.dim1() == 0 ? 0 : anysub_prob_neutral[b]);

      sim_info.anysub_prob = anysub_prob[b];
      sim_info.anysub_cdf = anysub_cdf[b];
      simulate_discrete_time_branch(treenode_seq[branch_node->parent()->node_id()],
                                    treenode_seq[branch_node->node_id()],
                                    etime,
                                    sim_info,
                                    sim_opt.is_report_time,
                                    apn_neutral,
                                    cat_sim_time[b]);
    }
  }

  // clear memory for internal node sequences, which we won't be using:
  //
  /// \todo simulation storage can be optimized down to
  /// seq_len*n_leaves if storage becomes a bigger issue
  //
  for(unsigned n(0);n<n_nodes;++n){
    if(! tree.node(n)->is_leaf()) treenode_seq[n].clear();
  }

  if(sim_opt.is_report_time){
    const unsigned bsize(SITE_MODEL::base_size(rg.site_model()));
    const smlfloat scale_factor(1./static_cast<smlfloat>(sim_opt.size*bsize));

    for(unsigned b(0);b<n_branches;++b){
      for(unsigned c(0);c<n_cats;++c){
        cat_sim_time[b][c] *= scale_factor;
      }
    }

    report_sim_tree_times(cm,tree,cat_sim_time,log_os);
  }

#ifdef DIRECT_SITE_CONV
  // write sim data out to sd:
  //
  sim_site_data_init(tree,mdl.rate_gtor_model(),sd);

  const SITE_MODEL::index_t smodel(RATE_GTOR_MODEL::convert_to_site_model(mdl.rate_gtor_model()));
  const unsigned ambig_state(SITE_MODEL::ambig_state(smodel));

  site_code sc(n_leaves);

  static const char group_label[] = "sim_discrete_group_";
  for(unsigned i(0);i<n_groups;++i){
    ostringstream oss;
    oss << group_label << i;
    sd.group_label.push_back(oss.str());
  }

  // chop up present day sequences into site_data data format:
  for(unsigned i(0);i<sim_opt.size;++i){

    for(unsigned t(0);t<n_leaves;++t){
      sc.set_taxid(t,treenode_seq[tree.leaf_node(t)->node_id()][i]);
    }

    unsigned cat_no(UNASSIGNED_CAT_ID);
    if(sim_opt.assigned_cat_prob>0.){
      const prob_t r = random_uniform();
      if(r<sim_opt.assigned_cat_prob){
        cat_no=sim_info.cat_seq[i]+1;
      }
    }

    const unsigned group_no(group_seq[i]);

    sd.site_count[sc][count_data_code(group_no,cat_no)]++;
  }
#else
  // write sim data out to nuc_seq_data:
  //
  const RATE_GTOR_MODEL::index_t rgm(mdl.rate_gtor_model());

  sim_seq_data_init(tree,n_groups,n_cats,"sim_discrete",cm,nsd);

  group_range gr(group_seq);

  const SITE_MODEL::index_t smodel(RATE_GTOR_MODEL::convert_to_site_model(rgm));
  const unsigned bsize(SITE_MODEL::base_size(smodel));

  unsigned codon_pos[SITE_MODEL::MAX_BASE_SIZE];
  SITE_MODEL::codon_position(smodel,codon_pos);

  for(unsigned g(0);g<n_groups;++g){
    unsigned group_start,group_stop;
    if(! gr.get_gbounds(g,group_start,group_stop)) continue;

    const unsigned group_size(group_stop-group_start);
    align_dat_t& na(nsd.dat[g]);
    na.seq.init(n_leaves,group_size*bsize);
    na.info.resize(group_size*bsize);

    for(unsigned i(0);i<group_size;++i){

      // relies on assumption that sites are non-overlapping:
      for(unsigned t(0);t<n_leaves;++t){
        const std::vector<unsigned>& leaf_seq(treenode_seq[tree.leaf_node(t)->node_id()]);
        SITE_MODEL::decode_nuc(smodel,leaf_seq[group_start+i],na.seq[t]+i*bsize);
      }

      unsigned site_cat_no(UNASSIGNED_CAT_ID);
      if(sim_opt.assigned_cat_prob>0.){
        const prob_t r(random_uniform());
        if(r<sim_opt.assigned_cat_prob){ site_cat_no=sim_info.cat_seq[group_start+i]+1; }
      }

      for(unsigned j(0);j<bsize;++j){
        const unsigned bi((i*bsize)+j);

        align_col_t& nai(na.info[bi]);

        nai.class_no=site_cat_no;
        nai.codon_pos=codon_pos[j];
        nai.is_continuous = (bi!=0);
      }
    }
  }
#endif
}

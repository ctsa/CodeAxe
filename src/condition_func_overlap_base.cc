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

// $Id: condition_func_overlap_base.cc 1148 2008-03-11 00:49:49Z ctsa $

/// \file

#include "bi_tree.h"
#include "condition_func_overlap_base.h"
#include "lhood_site_cached.h"
#include "lhood_site_cached_optorder.h"
#include "site_data_fastlup_core_util.h"
#include "subs_ml_ptol.h"
#include "tree_probs.h"
#include "util/math/prob_util.h"
#include "util/general/die.h"

#include <cstring>

#include <set>



void
condition_func_overlap_base::
data_init_cfd(cf_data& cfdx,
              const unsigned condition_no,
              const bi_tree& tree,
              const site_data_fastlup_core& sdf){

  const unsigned n_leaves(tree.leaf_size());
  const unsigned n_states(state_size());

  // include +1 ambiguous state
  std::vector<unsigned> state_translation_map_write(n_states+1);
  for(unsigned i(0);i<n_states+1;++i){
    state_translation_map_write[i] = state_2_conditioning_state(i,condition_no);
  }
  const std::vector<unsigned>& state_translation_map(state_translation_map_write);

  // rebuild sdf in the reduced conditioning_state space:
  //
  std::set<site_code> conditioned_code_set;
  for(unsigned i(0);i<sdf.len;++i){
    site_code sc(n_leaves);
    for(unsigned j(0);j<n_leaves;++j){
      sc.set_taxid(j,state_translation_map[sdf.data_core[i].index[j]]);
    }
    conditioned_code_set.insert(sc);
  }

  cfdx.sdf_cf.init(conditioned_code_set.size(),n_leaves);

  std::set<site_code>::const_iterator i=conditioned_code_set.begin(),i_end=conditioned_code_set.end();

  for(unsigned count(0);i!=i_end;++i,++count){
    for(unsigned j(0);j<n_leaves;++j){
      cfdx.sdf_cf.data_core[count].index[j] = i->get_taxid(j);
    }
  }

  // find optimal sdf_cf org order for cached lhood calculation
  lhood_site_cached_optorder(tree,n_states,cfdx.sdf_cf);

  cfdx.sdf_cf.fix();

  // phew! sdf_cf is done... let's initialize everything else:
  cfdx.data_reduced_state_index_map.resize(sdf.len);
  cfdx.condition_prob.init(cfdx.sdf_cf.len);

  map_translated_index_superset_data(sdf,cfdx.sdf_cf,state_translation_map.begin(),cfdx.data_reduced_state_index_map.begin());
}


void
condition_func_overlap_base::
data_init(const bi_tree& tree,
          const site_data_fastlup_core& sdf){

  if(tree_ptr){ die("multiple data_init in cf!!"); }
  tree_ptr=&tree;

  data_init_cfd(cfd,0,tree,sdf);
  if(is_bidirectional()){
    data_init_cfd(cfd2,1,tree,sdf);
  }
}



void
condition_func_overlap_base::
update_condition_prob(const cf_data& cfdx) const {

  const bi_tree& tree(*tree_ptr);
  const unsigned n_states(state_size());
  const unsigned n_conditioning_states(conditioning_state_size());

  const tree_lhood_info ti(tree,cfdx.tprob_cf,n_states,n_conditioning_states);

  lhood_site_cached(ti,cfdx.sdf_cf,cfdx.condition_prob.ptr());
}



void
condition_func_overlap_base::
tprob_init_cfd(const cf_data& cfdx,
               const unsigned condition_no,
               const tree_probs& tprob) const {

  const bi_tree& tree(*tree_ptr);
  const unsigned n_leaves(tree.leaf_size());
  const unsigned n_states(state_size());
  const unsigned n_conditioning_states(conditioning_state_size());

  // setup hacked tree_probs, with leaf branches replaced by
  // condition_func leaf_leg
  //
  cfdx.tprob_cf=tprob;

  for(unsigned b(0);b<n_leaves;++b){
    const int branch_id(tree.leaf_node(b)->branch_id());

    // initialize leaves using the class conditioning func
    memset(cfdx.tprob_cf.branch_prob(branch_id),0,n_states*n_conditioning_states*sizeof(prob_t));
    for(unsigned s(0);s<n_conditioning_states;++s){
      get_leaf_leg(cfdx.tprob_cf.branch_prob(branch_id)+s*n_states,
                   tprob.branch_prob(branch_id),s,condition_no);
    }

    // must set this to false (if true), because lhood algorithm
    // assumes an ident prob matrix at zero time.  because we're
    // simultaneously allowing for multiple states (all states which
    // map to the conditioning state), the leaf matrix will not be
    // ident even at zero time:
    //
    cfdx.tprob_cf.is_branch_tprob_identity(branch_id) = false;
  }

  update_condition_prob(cfdx);
}



void
condition_func_overlap_base::
tprob_init(const tree_probs& tprob) const {

  if(! tree_ptr){ die("cfob.condition_site_prob(): data not initialized"); }

  const bi_tree& tree(*tree_ptr);
  if(tree.branch_size() != tprob.n_branches) { die("cfob.condition_site_prob(): tree/tprob size mismatch"); }

#ifdef DEBUG
  tprob.self_check();
#endif

  tprob_init_cfd(cfd,0,tprob);
  if(is_bidirectional()){
    tprob_init_cfd(cfd2,1,tprob);
  }
}



void
condition_func_overlap_base::
root_update_cfd(const cf_data& cfdx,
                const prob_t* root) const {

  cfdx.tprob_cf.update_root(root);
  update_condition_prob(cfdx);
}



void
condition_func_overlap_base::
root_update(const prob_t* root) const {

  root_update_cfd(cfd,root);
  if(is_bidirectional()) root_update_cfd(cfd2,root);
}



void
condition_func_overlap_base::
condition_site_prob_prep(const std::vector<bool>& site_prob_mask,
                         prob_t* site_prob) const {

  // condition probs are calculated, now apply them to site prob:
  //
  const unsigned site_prob_size(cfd.data_reduced_state_index_map.size());

  for(unsigned i(0);i<site_prob_size;++i){
    if(site_prob_mask[i]) continue;
    smlfloat denom(cfd.condition_prob[cfd.data_reduced_state_index_map[i]]);
    if(is_bidirectional()){
      denom = std::sqrt(denom*cfd2.condition_prob[cfd2.data_reduced_state_index_map[i]]);
    }
    site_prob[i] /= denom;
  }
}

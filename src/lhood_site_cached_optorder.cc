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

// $Id: lhood_site_cached_optorder.cc 999 2007-11-02 18:48:00Z ctsa $

/// \file

#include "lhood_site_cached_optorder.h"
#include "site_code.h"
#include "util/general/die.h"

#include <algorithm>

#include <sstream>

#ifdef DEBUG
#include "util/general/log.h"

#include <ostream>
#endif

struct tree_lhood_info_optorder {

  const bi_tree& tree;
  site_data_fastlup_core& sdf;
  const int n_states;
};


static
void
lhood_org_order_score_combine_branches(const tree_lhood_info_optorder& ti,
                                       const bi_tree_node* current_node,
                                       std::vector<const bi_tree_node*>& prior_branch_lhood,
                                       const bi_tree_node*& branch_lhood_node,
                                       const bi_tree_node*& combined_lhood_node,
                                       unsigned& cost,
                                       int legal_next_leaf_id[2]){

  const bi_tree_node* next_node(0);

  // root case:
  if(current_node->is_root()){
    next_node=branch_lhood_node->sister();

  // test for prior case:
  } else {

    const bi_tree_node* edge_node[3] = { current_node->parent(),
                                         current_node->child1(),
                                         current_node->child2() };

    // put 2 continuing edges first in edge_node list
    for(unsigned i(0);i<2;++i){
      if(edge_node[i] == branch_lhood_node) {
        std::swap(edge_node[i],edge_node[2]);
      }
    }

    bool is_prior(false);
    for(unsigned i(0);i<prior_branch_lhood.size();++i){
      for(unsigned j(0);j<2;++j) {
        if(edge_node[j] == prior_branch_lhood[i]){
          is_prior=true;
          next_node=edge_node[(j+1)%2];
          break;
        }
      }
      if(is_prior) {
        std::swap(prior_branch_lhood[i],prior_branch_lhood.back());
        prior_branch_lhood.pop_back();
        break;
      }
    }

    // null case:
    if(! is_prior) {
      unsigned n_next(0);
      for(unsigned i(0);i<3;++i){
        if(edge_node[i] != branch_lhood_node){

          // first look for leaf nodes connected to current (up to 2)
          if(edge_node[i]->is_leaf()){
            legal_next_leaf_id[n_next] = edge_node[i]->leaf_id();
            n_next++;

            // also look for leaf node connected to current through the root:
          } else if(edge_node[i]->is_root()){
            const bi_tree_node* sister_node(current_node->sister());
            if(sister_node->is_leaf()){
              legal_next_leaf_id[n_next] = sister_node->leaf_id();
              n_next++;
            }
          }
        }
      }
      assert(n_next<3);
      return;
    }
  }

  combined_lhood_node=current_node;

  cost += ti.n_states;

  // check that we're not moving into a leaf node -- else time to
  // return node_cache back down the stack and move to next organism:
  if(! next_node->is_leaf()) {
    cost += ti.n_states*ti.n_states;
    branch_lhood_node=combined_lhood_node;
    lhood_org_order_score_combine_branches(ti,next_node,prior_branch_lhood,branch_lhood_node,combined_lhood_node,cost,legal_next_leaf_id);
  } else {
    // detect legal next leaf_id:
    legal_next_leaf_id[0] = next_node->leaf_id();
  }
}



struct minfo {
  minfo(const unsigned n_leaves) : is_found(false),cost(0),order(n_leaves) {}

  bool is_found;
  unsigned cost;
  std::vector<unsigned> order;
};


struct minfo2 {
  minfo2(const unsigned n_leaves) :
    min_cachemove_spectrum(n_leaves,n_leaves),
#if 0
    min_final_cost_upper(0), is_min_final_cost_upper(false),
#endif
    n_search(0) {}

  std::vector<unsigned> min_cachemove_spectrum;
#if 0
  unsigned min_final_cost_upper;
  bool is_min_final_cost_upper;
#endif
  unsigned n_search;
};


struct minfo3 {
  unsigned max_search;
  unsigned max_cachemove;
};


// recursive scoring function
//  \todo -- doc me --
static
void
lhood_org_order_score(tree_lhood_info_optorder& ti,
                      const std::vector<const bi_tree_node*>& prior_branch_lhood,
                      const bi_tree_node* prior_combined_lhood_node,
                      const int leaf_index,
                      unsigned cost,
                      minfo& mi,
                      minfo2& mi2,
                      const minfo3& mi3,
                      const unsigned remaining_branch_cost,
                      const unsigned n_cachemove = 0,
                      const bool is_last_move_cachemove = false){


  const unsigned leaf_id(ti.sdf.order[leaf_index]);
  const bi_tree_node* leaf_node(ti.tree.leaf_node(leaf_id));
  const bi_tree_node* parent_node(leaf_node->parent());

  if(parent_node==0) die("Can't handle singleton tree");

#ifdef DEBUG
  log_os << "entering lhood_oos: leaf_index: " << leaf_index << "\n";
  for(unsigned j(0);j<ti.sdf.order.size();++j){
    log_os << "sdf_order: " << j << " " << ti.sdf.order[j] << "\n";
  }
  //  log_os << "cost: " << cost << "\n";
#endif

  if(leaf_index+1 == static_cast<int>(ti.tree.leaf_size())){

    if(parent_node==prior_combined_lhood_node){
      // legal ordering:
      if(cost < mi.cost || ! mi.is_found){
        mi.cost = cost;
        mi.order = ti.sdf.order;
        mi.is_found = true;
#ifdef DEBUG
        for(unsigned j(0);j<mi.order.size();++j){
          log_os << "min_order: " << j << " " << mi.order[j] << "\n";
        }
        log_os << "cost: " << cost << "\n";
        log_os << "nnd: " << n_cachemove << "\n";
#endif
      }
    }
    return;
  }

  mi2.n_search += 1;
  if(mi3.max_search && mi.is_found && mi2.n_search > mi3.max_search) return;

  const bi_tree_node* branch_lhood_node(leaf_node);
  const bi_tree_node* combined_lhood_node(0);

  unsigned branch_cost(0);

  int legal_next_leaf_id[2] = {-1,-1};

  std::vector<const bi_tree_node*> local_prior_branch_lhood(prior_branch_lhood);

  // walk through tree to find cost of calculating on this branch
  lhood_org_order_score_combine_branches(ti,parent_node,local_prior_branch_lhood,
                                         branch_lhood_node,combined_lhood_node,
                                         branch_cost,legal_next_leaf_id);


  unsigned branch_count(0);
  //  if(branch_cost){
    std::vector<unsigned> leaf_state(leaf_index+1);
    for(unsigned i(0);i<ti.sdf.len;++i){
      bool is_count(false);
      for(int j(0);j<leaf_index+1;++j){
        const unsigned j_leaf_id(ti.sdf.order[j]);
        const unsigned new_leaf_state(ti.sdf.data_core[i].index[j_leaf_id]);
        if(i==0 || leaf_state[j] != new_leaf_state){
          leaf_state[j] = new_leaf_state;
          is_count=true;
        }
      }
      if(is_count) branch_count++;
    }
    //  }

  cost += branch_cost*branch_count;

  // solve for a lower and upper-bound of the final cost:
  const unsigned final_cost_lower(cost+branch_count*(remaining_branch_cost-branch_cost));
#if 0
  const unsigned final_cost_upper(cost+ti.sdf.len*(remaining_branch_cost-branch_cost));
#endif

  //  log_os << "leaf_index,branch_cost,branch_count,cost,final_cost: " << leaf_index << " " << branch_cost << " " << branch_count << " " << cost << " " << final_cost_lower << " " << final_cost_upper << "\n";

  // exit early if we've already exceeded the minimum cost:
  if(mi.is_found){
    if(final_cost_lower>=mi.cost) return;
  }
#if 0
  if(! mi2.is_min_final_cost_upper || final_cost_upper<mi2.min_final_cost_upper){
    mi2.min_final_cost_upper=final_cost_upper;
    mi2.is_min_final_cost_upper=true;
  }
  if(final_cost_lower>mi2.min_final_cost_upper) return;
#endif

  local_prior_branch_lhood.push_back(branch_lhood_node);

  const unsigned next_leaf_index(leaf_index+1);

  if(legal_next_leaf_id[0] != -1) { // regular move

    mi2.min_cachemove_spectrum[leaf_index] = n_cachemove;

    for(unsigned i(0);i<2;++i){
      const int swap_leaf_id(legal_next_leaf_id[i]);
      if(swap_leaf_id != -1){
        std::vector<unsigned>::iterator swap_leaf_id_it(std::find(ti.sdf.order.begin()+next_leaf_index,
                                                                  ti.sdf.order.end(),
                                                                  static_cast<unsigned>(swap_leaf_id)));
        if(swap_leaf_id_it == ti.sdf.order.end()){
          std::ostringstream oss;
          oss << __FILE__ << ":" << __LINE__ << " can't find swap_leaf_id index";
          throw substk_exception(oss.str().c_str());
        }
        const unsigned swap_leaf_index(swap_leaf_id_it-ti.sdf.order.begin());

        std::swap(ti.sdf.order[next_leaf_index],ti.sdf.order[swap_leaf_index]);
        if(next_leaf_index != swap_leaf_index || i!=0) {
          std::sort(ti.sdf.data_core,ti.sdf.data_core+ti.sdf.len,ti.sdf.cell_sorter);
        }

        lhood_org_order_score(ti,local_prior_branch_lhood,combined_lhood_node,next_leaf_index,cost,mi,mi2,mi3,
                              remaining_branch_cost-branch_cost,n_cachemove,false);

        std::swap(ti.sdf.order[next_leaf_index],ti.sdf.order[swap_leaf_index]);
      }
    }
  } else {  // cache move
    //exit early for singleton, no-legal next branch:
    if(branch_cost==0) return;

    // no two cachemoves in a row:
    if(is_last_move_cachemove) return;

    if(n_cachemove >= mi3.max_cachemove) return;

    // already found a solution which has fewer cached moves...
    if(n_cachemove >= mi2.min_cachemove_spectrum[leaf_index]) return;

    mi2.min_cachemove_spectrum[leaf_index] = n_cachemove+1;

    // no legal next leaf... start on next leaf to pick up the branch cache list elsewhere:
    //
    for(unsigned i(next_leaf_index);i<ti.sdf.order.size();++i){
      std::swap(ti.sdf.order[next_leaf_index],ti.sdf.order[i]);
      if(next_leaf_index != i) {
        std::sort(ti.sdf.data_core,ti.sdf.data_core+ti.sdf.len,ti.sdf.cell_sorter);
      }

      lhood_org_order_score(ti,local_prior_branch_lhood,combined_lhood_node,next_leaf_index,cost,mi,mi2,mi3,
                            remaining_branch_cost-branch_cost,n_cachemove+1,true);

      std::swap(ti.sdf.order[next_leaf_index],ti.sdf.order[i]);
    }
  }
}




// entry point to score function:
// scores the complexity of leaf order given a tree and site set
static
void
lhood_org_order_score(tree_lhood_info_optorder& ti,
                      minfo& mi,
                      minfo2& mi2,
                      const minfo3& mi3,
                      const unsigned total_branch_cost){

  static const unsigned cost(0);
  static const int leaf_index(0);
  std::vector<const bi_tree_node*> prior_branch_lhood;

  lhood_org_order_score(ti,prior_branch_lhood,0,leaf_index,cost,mi,mi2,mi3,total_branch_cost);
}



void
lhood_site_cached_optorder(const bi_tree& tree,
                           const int n_states,
                           site_data_fastlup_core& sdf){

  const unsigned n_leaves(tree.leaf_size());
  const unsigned n_nodes(tree.node_size());
  const unsigned n_branches(tree.branch_size());
  const unsigned n_nl_nodes(n_nodes-n_leaves);
  const unsigned n_nl_branches(n_branches-n_leaves);
  const unsigned total_branch_cost(n_states*(n_nl_nodes+n_states*n_nl_branches));

  minfo mi(n_leaves);
  minfo2 mi2(n_leaves);
  minfo3 mi3;

  tree_lhood_info_optorder ti = {tree,sdf,n_states};

  // start with naive sdf.order:
  for(unsigned j(0);j<n_leaves;++j){ sdf.order[j] =j; }

  // iterate through a few values of max_cachemove, then give up
  // and take a suboptimal tree if an ordering is not yet found:
  // (yes, this is dumb -- might be possible to rephrase this as
  // TSP and apply a standard near optimal solver to that --
  // see Kosakovsky-Pond & Muse)
  //
  mi3.max_cachemove = 0;
  mi3.max_search = 0;
  while(true){

    // iterate through starting leaves:
    for(unsigned i(0);i<n_leaves;++i){
#ifdef DEBUG
      log_os << "max_cachemove: " << mi3.max_cachemove << "\n";
      log_os << "sdforder, starting leaf: " << i << "\n";
#endif
      // swap in first position:
      std::swap(sdf.order[0],sdf.order[i]);
      std::sort(sdf.data_core,sdf.data_core+sdf.len,sdf.cell_sorter);

      lhood_org_order_score(ti,mi,mi2,mi3,total_branch_cost);

      // restore order:
      std::swap(sdf.order[0],sdf.order[i]);
    }

    if(mi.is_found){
      break;
    } else {
      // \todo -- doc me --
      //
      if(mi3.max_cachemove>=4){
        if(mi3.max_cachemove>=n_leaves){
          die("No legal org_ordering found!");
        }
        // give up and take a suboptimal tree
        mi3.max_cachemove=n_leaves;
        mi3.max_search=10000;
        mi2.n_search=0;
      } else {
        mi3.max_cachemove += 1;
      }
    }
  }

  sdf.order=mi.order;

#ifdef DEBUG
  for(unsigned j(0);j<n_leaves;++j){
    log_os << "min_order: " << j << " " << sdf.order[j] << "\n";
  }
#endif
}

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

// $Id: sim_util.cc 1193 2008-03-29 03:14:54Z ctsa $

/// \file

#include "bi_tree.h"
#include "cat_info.h"
#include "cat_manager.h"
#include "rate_gtor.h"
#include "report_util.h"
#include "root_gtor.h"
#include "sim_options.h"
#include "sim_util.h"
#include "subs_ml_model.h"
#include "simple_util.h"
#include "util/math/prob_util.h"
#include "util/math/random_util.h"

#include <sstream>


const smlfloat SIMULATE_DISCRETE_TIME_UNIT(5e-5);


void
get_ancestral_seq(const subs_ml_model& mdl,
                  const unsigned size,
                  const std::vector<unsigned>& cat_seq,
                  std::vector<unsigned>& seq){

  const unsigned n_states(mdl.state_size());
  const unsigned n_cats(mdl.get_cat_manager().cat_size());

  simple_matrix<prob_t> cdf(n_cats,n_states);
  for(unsigned i(0);i<n_cats;++i){
    pdistro_to_cdf(mdl.get_root_gtor().cat_state_pdistro(i),cdf[i],n_states);
  }

  seq.resize(size);
  for(unsigned i(0);i<size;++i) {
    seq[i] = random_cdf_variate(cdf[cat_seq[i]],n_states);
  }
}



/// 1. assign sites to groups:
///
void
get_sim_group_seq(std::vector<unsigned>& group_seq,
                  const sim_options& sim_opt){

  group_seq.resize(sim_opt.size);

  const unsigned group_range(sim_opt.group_size_max-sim_opt.group_size_min);
  unsigned group_no(0);
  unsigned current(0);
  while(current<sim_opt.size){
    unsigned group_len(sim_opt.group_size_min+static_cast<unsigned>(random_uniform()*static_cast<double>(group_range)));
    const int free_size(sim_opt.size-(current+group_len));
    if(free_size<static_cast<int>(sim_opt.group_size_min)) group_len = sim_opt.size-current;
    const unsigned group_end(current+group_len);
    for(;current<group_end;++current) group_seq[current] = group_no;

    group_no++;
  }
}



void
get_sim_cat_seq(std::vector<unsigned>& cat_seq,
                const std::vector<unsigned>& group_seq,
                const unsigned size,
                const subs_ml_model& mdl){

  cat_seq.resize(size);

  const cat_manager& cm(mdl.get_cat_manager());

  const unsigned n_cats(cm.cat_size());
  const unsigned n_group_cats(cm.group_cat_size());

  simple_array<prob_t> cat_pdistro(n_cats);
  cm.cat_pdistro(cat_pdistro.ptr());
  pdistro_to_cdf_inplace(cat_pdistro.ptr(),n_cats);

  simple_array<prob_t> group_cat_pdistro(n_group_cats);
  cm.group_cat_pdistro(group_cat_pdistro.ptr());
  pdistro_to_cdf_inplace(group_cat_pdistro.ptr(),n_group_cats);

  std::vector<unsigned> cat_group_switch(n_group_cats+1);
  {
    unsigned group_cat(0);
    for(unsigned i(0);i<n_cats;++i){
      if(cm.group_cat_no(i)==group_cat){
        cat_group_switch[group_cat]=i;
        group_cat++;
      }
    }
    assert(group_cat==n_group_cats);
    cat_group_switch[n_group_cats]=n_cats;
  }


  unsigned current_group_cat(0);
  for(unsigned i(0);i<size;++i){
    if(i==0 || group_seq[i-1] != group_seq[i]){
      current_group_cat = random_cdf_variate(group_cat_pdistro.ptr(),n_group_cats);
    }
    cat_seq[i] = random_subset_cdf_variate(cat_pdistro.ptr(),
                                           cat_group_switch[current_group_cat],
                                           cat_group_switch[current_group_cat+1]);

    assert(cm.group_cat_no(cat_seq[i])==current_group_cat);
  }
}


#if 0
void
sim_site_data_init(const bi_tree& tree,
                   const RATE_GTOR_MODEL::index_t rgm,
                   site_data& sd){

  sd.clear();

  const unsigned n_leaves(tree.leaf_size());

  for(unsigned i(0);i<n_leaves;++i){
    sd.taxid.assignid(tree.leaf_node(i)->label());
  }
  sd.sm = RATE_GTOR_MODEL::convert_to_site_model(rgm);

  /// \todo would prefer that we just write "simi" for i=1..N classes
  sd.data_class_labels.assignid(UNASSIGNED_CAT_LABEL);
}
#endif



void
sim_seq_data_init(const bi_tree& tree,
                  const unsigned n_groups,
                  const unsigned n_cats,
                  const char* const gtag,
                  const cat_manager& cm,
                  nuc_seq_data& nsd){

  nsd.clear();

  const unsigned n_leaves(tree.leaf_size());

  for(unsigned i(0);i<n_leaves;++i){
    nsd.taxid.assignid(tree.leaf_node(i)->label());
  }

  nsd.data_class_labels.assignid(std::string(gtag)+"_seq_cat_mixture");
  for(unsigned i(0);i<n_cats;++i){
    std::ostringstream oss;
    oss << gtag << "_seq_cat_" << cm.cat_label(i);
    nsd.data_class_labels.assignid(oss.str());
  }

  static const char group_label[] = "_group_";
  nsd.dat.resize(n_groups);
  for(unsigned i(0);i<n_groups;++i){
    std::ostringstream oss;
    oss << gtag << group_label << i;
    nsd.dat[i].label=oss.str();
  }
}



bool
group_range::
get_gbounds(const unsigned group_no,
            unsigned& gstart,
            unsigned& gstop){

  for(gstart=_lgs;_gs[gstart]<group_no;gstart++){};

  if(_gs[gstart] != group_no) return false;

  for(gstop=gstart;(gstop!=_gs.size() && _gs[gstop]<=group_no);gstop++){};

  gstart += _lp;
  if(gstop<_rp) return false;
  gstop -= _rp;

  if(gstart>=gstop) return false;
  _lgs=gstop;

  return true;
}



void
report_sim_tree_times(const cat_manager& cm,
                      const bi_tree& tree,
                      const std::vector<std::vector<smlfloat> >& cat_sim_time,
                      std::ostream& os){

  const unsigned n_cats(cm.cat_size());
  const unsigned n_tree_cats(cm.typed_cat_size(CAT_PARAM_TYPE::TIME));
  const unsigned n_branches(cm.branch_size());

  std::vector<smlfloat> expected_branch_time(n_branches);
  std::vector<std::vector<smlfloat> > tree_cat_branch_time(n_tree_cats);

  for(unsigned i(0);i<n_tree_cats;++i){
    tree_cat_branch_time[i].resize(n_branches);
  }

  for(unsigned i(0);i<n_branches;++i){
    expected_branch_time[i] = 0.;

    for(unsigned j(0);j<n_tree_cats;++j){
      tree_cat_branch_time[j][i] = 0.;
    }

    for(unsigned j(0);j<n_cats;++j){
      expected_branch_time[i] += cat_sim_time[i][j];

      const unsigned tree_cat_no(cm.typed_cat_no_from_cat_no(j,CAT_PARAM_TYPE::TIME));
      tree_cat_branch_time[tree_cat_no][i] += cat_sim_time[i][j];
    }
  }

  simple_array<prob_t> tree_cat_pdistro(n_tree_cats);
  cm.typed_cat_pdistro(tree_cat_pdistro.ptr(),CAT_PARAM_TYPE::TIME);
  for(unsigned i(0);i<n_tree_cats;++i){
    for(unsigned j(0);j<n_branches;++j){
      tree_cat_branch_time[i][j] /= tree_cat_pdistro[i];
    }
  }

  os << "\n";
  report_time_instance("expected_sim_time",expected_branch_time,tree,os);
  if(n_tree_cats>1){
    for(unsigned i(0);i<n_tree_cats;++i){
      std::ostringstream oss;
      oss << "tree_cat_" << i << "_sim_time";
      report_time_instance(oss.str().c_str(),tree_cat_branch_time[i],tree,os);
    }
  }
}

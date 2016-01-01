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

// $Id: site_data_fastlup_util.cc 1182 2008-03-27 01:55:02Z ctsa $

/// \file

#include "bi_tree.h"
#include "cat_info.h"
#include "cat_manager.h"
#include "data_class_assignment_util.h"
#include "lhood_site_cached_optorder.h"
#include "site_data.h"
#include "site_data_fastlup_util.h"
#include "util/general/log.h"

#include <algorithm>
#include <ostream>
#include <set>


// transfer info from map to fixed array
//
static
void
site_data_to_site_data_fastlup(const bi_tree& tree,
                               const unsigned n_states,
                               const data_class_id_assignment_map& dc_ads_map,
                               const site_data& sd,
                               site_data_fastlup& sdf){

  const unsigned n_seqs(tree.leaf_size());

  if(sd.taxid.size() < n_seqs){
    log_os << "ERROR:: fewer seqs in data ("<< sd.taxid.size() << ") than in tree (" << n_seqs << ")\n";
    abort();
  }


  // group_state_reduction (renumber group ids into smallest
  // contiguous set from 0):
  std::map<site_data_code::group_id_type,site_data_fastlup_code::group_id_type> group_sd_to_sdf_map;

  {
    sdf.group_label.clear();
    site_data_fastlup_code::group_id_type group_id(0);

    site_count_map::const_iterator i,i_end=sd.site_count.end();
    for(i=sd.site_count.begin();i!=i_end;++i){
      site_count_data_type::const_iterator j,j_end=i->second.end();
      for(j=i->second.begin();j!=j_end;++j){
        // mapping is 1 to 1:
        const site_data_code::group_id_type& gi(j->first.group_id);
        if(group_sd_to_sdf_map.count(gi) == 0){
          group_sd_to_sdf_map[gi] = group_id++;

          // transfer group_labels into reorderd set:
          sdf.group_label.push_back(sd.group_label[gi]);
        }
      }
    }
  }

  const unsigned n_sites(sd.site_count.size());
  const unsigned n_adsets(dc_ads_map.assigned_data_set_size());
  sdf.init(n_sites,n_seqs,n_adsets);

  // used to track which skiped data_classes we've already warned on:
  std::set<unsigned> is_data_class_warning;

  std::vector<unsigned> sd_org_tree_order(n_seqs);
  for(unsigned j(0);j<n_seqs;++j){
    sd_org_tree_order[j] = sd.taxid.getid(tree.leaf_node(j)->label());
  }

  unsigned n(0);
  site_count_map::const_iterator i,i_end(sd.site_count.end());
  for(i=sd.site_count.begin();i!=i_end;++i,++n){
    const site_code& sc(i->first);
    const site_count_data_type& site_count_data(i->second);

    // convenient ref to current sdf cell:
    site_data_fastlup_cell& dn(sdf.data[n]);

    // transfer site type to current sdf cell:
    for(unsigned j(0);j<n_seqs;++j){
      dn.index[j] = sc.get_taxid(sd_org_tree_order[j]);
    }

    std::map<site_data_fastlup_code,unsigned> group_assigned_data_set_map;

    // iterate through all group/data_class counts for this site type:
    site_count_data_type::const_iterator j,j_end(site_count_data.end());
    for(j=site_count_data.begin();j!=j_end;++j){
      const site_data_fastlup_code::group_id_type& group_id(group_sd_to_sdf_map[j->first.group_id]);
      const site_data_code::data_class_id_type& data_class_id(j->first.data_class_id);
      const unsigned& count(j->second);

      if(! dc_ads_map.is_data_class_assigned(data_class_id)){
        if(is_data_class_warning.count(data_class_id) == 0){
          const std::string& dl(sd.data_class_labels.getstr(data_class_id));
          log_os << "WARNING:: discarding unmapped data class: " << dl << "\n";
          is_data_class_warning.insert(data_class_id);
        }
        continue;
      }

      const unsigned assigned_data_set_no(dc_ads_map.get_data_class_assigned_data_set_id(data_class_id));
      group_assigned_data_set_map[site_data_fastlup_code(group_id,assigned_data_set_no)] += count;
    }

    dn.group_assigned_data_set_count.insert(dn.group_assigned_data_set_count.end(),
                                            group_assigned_data_set_map.begin(),
                                            group_assigned_data_set_map.end());
  }

  // find optimal sdf org order for cached lhood calculation:
  lhood_site_cached_optorder(tree,n_states,sdf);

  sdf.fix();

#if 0
  for(unsigned j(0);j<n_sites;++j){
    log_os << "site: " << j;
    for(unsigned s(0);s<n_seqs;++s){
      log_os << " " << sdf.data_core[j].index[s];
    }
    log_os << "\n";
  }
#endif
}



void
mdl_site_data_fastlup_init(const site_data& sd,
                           const subs_ml_model& mdl,
                           site_data_fastlup& sdf){

  data_class_id_assignment_map dc_ads_map(mdl.get_cat_manager(),sd);

  site_data_to_site_data_fastlup(mdl.tree(),
                                 mdl.state_size(),
                                 dc_ads_map,
                                 sd,sdf);
}



void
site_data_fastlup_state_reduction(const bi_tree& tree,
                                  const unsigned reduced_n_states,
                                  const site_data_fastlup& sdf_in,
                                  const unsigned* index_translation,
                                  site_data_fastlup& sdf_out){

  /// \todo fix this so that we don't have to loop back through
  /// site_data struct
  ///

  const unsigned n_adsets(sdf_in.n_assigned_data_sets);
  const unsigned n_seqs(sdf_in.n_orgs);

  site_data sd;

  // insert remapped counts into sd
  for(unsigned i(0);i<sdf_in.len;++i){
    site_code sc(n_seqs);
    for(unsigned t(0);t<n_seqs;++t){
      sc.set_taxid(t,index_translation[sdf_in.data[i].index[t]]);
    }

    const unsigned asize(sdf_in.data[i].group_assigned_data_set_count.size());
    for(unsigned j(0);j<asize;++j){
      const unsigned group_id(sdf_in.data[i].group_assigned_data_set_count[j].first.group_id);
      const unsigned adset_id(sdf_in.data[i].group_assigned_data_set_count[j].first.assigned_data_set_id);
      const unsigned count(sdf_in.data[i].group_assigned_data_set_count[j].second);

      // adset_id is being used here as the data class for the hack site_data structure:
      sd.site_count[sc][site_data_code(group_id,adset_id)] += count;
    }
  }

  // sloppy hack -- naive data_class_id_assignment:
  data_class_id_assignment_map dc_ads_map(n_adsets);

  // give sd taxid's that match the tree order:
  std::vector<unsigned> sd_org_tree_order(n_seqs);
  for(unsigned j(0);j<n_seqs;++j){
    sd.taxid.assignid(tree.leaf_node(j)->label());
  }

  site_data_to_site_data_fastlup(tree,reduced_n_states,dc_ads_map,sd,sdf_out);
}

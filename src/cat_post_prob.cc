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

// $Id: cat_post_prob.cc 1182 2008-03-27 01:55:02Z ctsa $

/// \file

#include "bi_tree.h"
#include "cat_manager.h"
#include "cat_post_prob.h"

#include <cmath>

#include <ostream>
#include <string>
#include <vector>




void
make_cat_post_prob(const unsigned n_cats,
                   const site_data_fastlup& full_sdf,
                   const CPPM::index_t cat_post_prob_mode,
                   prob_t** ppi){

  if(full_sdf.n_assigned_data_sets > 1){
    die("cat posterior prob functions not yet enabled for multiple data sets");
  }

  if       (cat_post_prob_mode == CPPM::SITE){
    for(unsigned i(0);i<full_sdf.len;++i){
      smlfloat sum(0.);
      for(unsigned c(0);c<n_cats;++c) sum += ppi[i][c];
      for(unsigned c(0);c<n_cats;++c) ppi[i][c] /= sum;
    }
  } else if(cat_post_prob_mode == CPPM::GROUP){
    const unsigned n_groups(full_sdf.n_groups);

    for(unsigned i(0);i<n_groups;++i){
      const smlfloat mscale(*(std::max_element(ppi[i],ppi[i]+n_cats)));
      for(unsigned c(0);c<n_cats;++c){
        ppi[i][c] = std::exp(ppi[i][c]-mscale);
      }

      smlfloat sum(0.);
      for(unsigned c(0);c<n_cats;++c) sum += ppi[i][c];
      for(unsigned c(0);c<n_cats;++c) ppi[i][c] /= sum;
    }
  }
}



void
cat_post_prob_report(const RATE_GTOR_MODEL::index_t rgm,
                     const bi_tree& tree,
                     const cat_manager& cm,
                     const site_data_fastlup& full_sdf,
                     const prob_t * const * ppi,
                     const CPPM::index_t cat_post_prob_mode,
                     std::ostream& os){

  assert(full_sdf.n_orgs == tree.leaf_size());

  const SITE_MODEL::index_t sm(RATE_GTOR_MODEL::convert_to_site_model(rgm));

  const unsigned n_cats(cm.cat_size());

  // first write out site info:
  os << "base_size: " << SITE_MODEL::base_size(sm) << "\n";
  os << "repeat_size: " << RATE_GTOR_MODEL::base_size_conditioned(rgm) << "\n";
  os << "repeat_offset: " << RATE_GTOR_MODEL::base_size_repeat_offset(rgm) << "\n";

  os << "\n";
  std::vector<prob_t> cat_pdistro(n_cats);
  cm.cat_pdistro(cat_pdistro.begin());
  for(unsigned c(0);c<n_cats;++c){
    os << "prior_prob: " << cm.cat_label(c) << " " << cat_pdistro[c] << "\n";
  }
  os << "\n";

  if       (cat_post_prob_mode == CPPM::SITE){
    const unsigned n_leaves(tree.leaf_size());

    // write out org id's:
    for(unsigned j(0);j<n_leaves;++j){
      if(j) os << " ";
      os << tree.leaf_node(j)->label();
    }

    for(unsigned c(0);c<n_cats;++c){
      os << " " << cm.cat_label(c);
    }
    os << "\n";

    for(unsigned i(0);i<full_sdf.len;++i){
      for(unsigned j(0);j<n_leaves;++j){
        if(j) os << " ";
        os << SITE_MODEL::print(sm,full_sdf.data_core[i].index[j]);
      }

      for(unsigned c(0);c<n_cats;++c) os << " " << ppi[i][c];
      os << "\n";
    }
  } else if(cat_post_prob_mode == CPPM::GROUP){
    const unsigned n_groups(full_sdf.n_groups);

    os << "group_label";
    for(unsigned c(0);c<n_cats;++c){
      os << " " << cm.cat_label(c);
    }
    os << "\n";

    for(unsigned i(0);i<n_groups;++i){
      os << full_sdf.group_label[i];
      for(unsigned c(0);c<n_cats;++c) os << " " << ppi[i][c];
      os << "\n";
    }
  }
}

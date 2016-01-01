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

// $Id: report_util.hh 918 2007-10-12 23:36:13Z ctsa $

/// \file

#include "util/math/indy_random_var_io.h"

#include <ostream>


template <typename FloatType>
void
report_time_instance(const char* time_label,
                     const std::vector<FloatType>& branch_time,
                     const bi_tree& tree,
                     std::ostream& os) {

  const unsigned n_branches(tree.branch_size());

  for(unsigned i(0);i<n_branches;++i){
    os << time_label << ": " << tree.branch_node(i)->label() << " : " << branch_time[i] << "\n";
  }
  os << "\n";

  //also report time as a newick tree
  bi_tree t2(tree);
  for(unsigned i(0);i<n_branches;++i){
    t2.branch_node(i)->set_length(static_cast<smlfloat>(branch_time[i]));
  }
  os << t2 << "\n\n";
}

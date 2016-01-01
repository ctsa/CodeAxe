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

// $Id: process_fsa.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __PROCESS_FSA_H
#define __PROCESS_FSA_H


#include "nuc_seq_data.h"


struct seq_filter_options {
  seq_filter_options() :
    is_use_window_filter(true),
    is_use_ambi_filter(true),
    is_use_gap_filter(true),
    is_skip_filters(false),
    min_match(0),
    max_align_count(0) {}

  bool is_use_window_filter;
  bool is_use_ambi_filter;
  bool is_use_gap_filter;
  bool is_skip_filters;
  unsigned min_match;
  unsigned max_align_count;
};


/// \brief get filtered nuc_seq_data from (dir of) seq file(s)
///
void
process_seq_data(nuc_seq_data& in_seq,
                 const seq_filter_options& fopt,
                 const char* seq_file_node,
                 const char* cat_file_node = 0,
                 const char* cat_labels_file = 0);

#endif

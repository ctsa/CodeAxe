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

// $Id: nuc_seq_data.h 1193 2008-03-29 03:14:54Z ctsa $

/// \file

#ifndef __NUC_SEQ_DATA_H
#define __NUC_SEQ_DATA_H


#include "audit_info.h"
#include "name_id_lup.h"
#include "simple_util.h"
#include "util/bio/bioseq_util.h"

#include <vector>


extern const unsigned DEFAULT_CLASS_ID;


struct align_col_t {

  typedef unsigned char class_type;

  align_col_t() : class_no(DEFAULT_CLASS_ID), codon_pos(0), is_continuous(true) {}

  unsigned char class_no;
  unsigned char codon_pos;
  bool is_continuous;       ///< is this column continuous with the preceding column?
};


struct align_dat_t {
  typedef std::vector<align_col_t> info_type;
  typedef simple_init_matrix<NUC::index_t> seq_type;

  std::string label;
  seq_type seq; // n_seqs*ncol
  info_type info;
};



/// \brief intermediate data format between seq file/dir read and
/// conversion to independent site data exchange format
///
/// this is an abstraction of sequence and data-class information from any file
/// format, and will be the default format for seq-dependent methods (MCMC)
///
struct nuc_seq_data {

  nuc_seq_data() : _ai() {}

  explicit
  nuc_seq_data(const audit_info& ai) : _ai(ai) {}

  unsigned n_seqs() const { return taxid.size(); }

  unsigned data_class_size() const { return data_class_labels.size(); }

  void load_state(std::istream& is);

  void store_state(std::ostream& os) const;

  void clear() {
    data_class_labels.clear();
    taxid.clear();
    dat.clear();
  }

  name_id_lup data_class_labels;
  name_id_lup taxid;
  std::vector<align_dat_t> dat;

private:
  const audit_info _ai;
};


void
nuc_seq_data_report(const nuc_seq_data& in_seq,
                    std::ostream& os);

#endif

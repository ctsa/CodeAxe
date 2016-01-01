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

// $Id: nuc_seq_data.cc 1193 2008-03-29 03:14:54Z ctsa $

/// \file

#include "nuc_seq_data.h"
#include "subs_ml_types.h"
#include "util/bio/bioseq_util.h"
#include "util/general/log.h"

#include <cassert>

#include <iostream>
#include <sstream>


const unsigned DATA_FILE_VERSION(1);

const unsigned DEFAULT_CLASS_ID(0);


void
nuc_seq_data_report(const nuc_seq_data& in_seq,
                    std::ostream& os){

  const unsigned n_alignments(in_seq.dat.size());
  for(unsigned i(0);i<n_alignments;++i){
    const align_dat_t::seq_type& seq(in_seq.dat[i].seq);
    const unsigned n_seqs(seq.dim1());
    const unsigned n_bases(seq.dim2());

    unsigned gc_count(0),total(0);
    for(unsigned t(0);t<n_seqs;++t){
      const NUC::index_t* n(seq[t]);
      for(unsigned j(0);j<n_bases;++j){
        if( n[j] == NUC::G || n[j] == NUC::C ) gc_count++;
        if( n[j] != NUC::N ) total++;
      }
    }

    if(total){
      const smlfloat gc(static_cast<smlfloat>(gc_count)/static_cast<smlfloat>(total));
      os << "GENE_GC: " << i << " " << total/n_seqs << " " << gc << "\n";
    }
  }
}

static const char* org_signal("org:");
static const char* class_signal("data_class:");
static const char align_signal('>');
static const char section_signal(']');



const char* const version_tag = "nuc_seq_file_version";


void
nuc_seq_data::
store_state(std::ostream& os) const {
  os << version_tag << " " << DATA_FILE_VERSION << "\n";
  _ai.store_state(os);

  const unsigned ns(n_seqs());

  os << org_signal;
  for(unsigned t(0);t<ns;++t){
    os << " " << taxid.getstr(t);
  }
  os << "\n";

  const unsigned n_classes(data_class_size());

  os << class_signal;
  for(unsigned i(0);i<n_classes;++i){
    os << " " << data_class_labels.getstr(i);
  }
  os << "\n";

  const unsigned n_alignments(dat.size());
  for(unsigned i(0);i<n_alignments;++i){
    const align_dat_t& d(dat[i]);
    const unsigned n_bases(d.seq.dim2());
    os << align_signal << " " << d.label << " " << n_bases << "\n";

    for(unsigned j(0);j<n_bases;++j){
      if(! d.info[j].is_continuous) {
        os << section_signal << static_cast<unsigned>(d.info[j].codon_pos) << "\n";
      }

      for(unsigned k(0);k<ns;++k){
        os << d.seq[k][j];
      }
      if(d.info[j].class_no != DEFAULT_CLASS_ID) os << static_cast<unsigned>(d.info[j].class_no);
      os << "\n";
    }
  }
}



static
void
load_die(std::istream&,
         const char* const sid,
         const std::string& tmpstr){

  log_os << "ERROR:: -- Invalid nuc sequence data file format. --\n"
         << "ERROR:: expected_section: " << sid << "\n"
         << "ERROR:: section: " << tmpstr << "\n"
         << "\n";
  //  log_os << is.rdbuf();
  abort();
}




static
void
load_advance(std::istream& is,
             const char* const second = 0){

  std::string tmpstr;
  is >> tmpstr;
  if(second && tmpstr != second) load_die(is,second,tmpstr);
}


void
nuc_seq_data::
load_state(std::istream& is) {

  clear();

  {  // check file version
    load_advance(is,version_tag);
    unsigned tmpu;
    is >> tmpu;
    if(tmpu != DATA_FILE_VERSION){
      log_os << "ERROR:: bad nuc sequence data file version number. expected: " << DATA_FILE_VERSION << "\n";
      exit(EXIT_FAILURE);
    }
  }

  _ai.skip_state(is);

  {
    std::string buf;

    while(is){
      if(is.peek()==align_signal) break;

      if( ! getline(is,buf) ) break;

      std::istringstream sbuf(buf);
      std::string tmpstr;

      sbuf >> tmpstr;
      if(tmpstr==org_signal){
        while(sbuf >> tmpstr){ taxid.assignid(tmpstr); }
      }
      if(tmpstr==class_signal){
        while(sbuf >> tmpstr) { data_class_labels.assignid(tmpstr); }
      }
    }
  }

  const unsigned ns(n_seqs());

  bool is_new_segment(false);
  unsigned codon_pos_sync(0);
  unsigned base_no(0);
  while(is){
    int c=is.get();
    if(c==' ' || c=='\t' || c== '\n' || ! is) continue;
    if(c==align_signal){

      dat.push_back(align_dat_t());

      unsigned n_bases;
      is >> dat.back().label >> n_bases;

      dat.back().info.resize(n_bases);
      dat.back().seq.init(ns,n_bases);

      base_no=0;

    } else if(c==section_signal){
      is_new_segment=true;

      if(is.peek() != '\n'){
        unsigned tmpu;
        is >> tmpu;
        codon_pos_sync=((tmpu+CODON::BASE_SIZE)-base_no%CODON::BASE_SIZE)%CODON::BASE_SIZE;
      } else {
        codon_pos_sync=0;
      }

    } else {
      for(unsigned i(0);i<ns;++i){
        if(i!=0) c=is.get();
        assert(c != ' ' && c != '\n');
        assert(is);
        dat.back().seq[i][base_no] = static_cast<NUC::index_t>(c-'0');
      }

      if(is.peek() != '\n'){
        unsigned tmpu;
        is >> tmpu;
        dat.back().info[base_no].class_no=static_cast<unsigned char>(tmpu);
      }

      dat.back().info[base_no].codon_pos =
        static_cast<unsigned char>((base_no+codon_pos_sync)%CODON::BASE_SIZE);
      dat.back().info[base_no].is_continuous = !is_new_segment;

      is_new_segment=false;
      base_no++;
    }
  }
}

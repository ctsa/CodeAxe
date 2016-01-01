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

// $Id: site_data_util.cc 1193 2008-03-29 03:14:54Z ctsa $

/// \file

#include "cat_info.h"
#include "nuc_seq_data.h"
#include "site_data.h"
#include "site_data_util.h"
#include "util/bio/bioseq_util.h"
#include "util/general/log.h"

#include <ostream>
#include <vector>

using namespace std;


static
void
make_gc_include_map(std::vector<bool>& include_seq,
                    const nuc_seq_data& in_seq,
                    const double gc_max,
                    const double gc_min){

  const unsigned n_alignments(in_seq.dat.size());
  include_seq.resize(n_alignments);
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

    include_seq[i] = false;
    if(total){
      const double gc(static_cast<double>(gc_count)/static_cast<double>(total));
      if( gc <= gc_max && gc > gc_min ) include_seq[i] = true;
    }
  }
}




static
void
check_codon_data_class(const align_dat_t& in_align,
                       const unsigned base_no){

  if( in_align.info[base_no].class_no != in_align.info[base_no+1].class_no ||
      in_align.info[base_no].class_no != in_align.info[base_no+2].class_no ){
    log_os << "ERROR:: data class changes within codon\n";
    log_os << "cat@pos1,2,3: "
           << in_align.info[base_no].class_no << " "
           << in_align.info[base_no+1].class_no << " "
           << in_align.info[base_no+2].class_no << "\n";
    abort();
  }
}


// take new site_code and add it to site_data, with cat info
static
void
store_site_info(site_data& sd,
                const align_dat_t& in_align,
                const unsigned align_no,
                const unsigned base_no,
                const site_code& sc,
                const bool is_codon_type){

  if(is_codon_type) {
    check_codon_data_class(in_align,base_no);
  }

  const unsigned data_class_id(in_align.info[base_no].class_no);

  sd.site_count[sc][site_data_code(align_no,data_class_id)]++;
}



static
void
process_site(site_data& sd,
             const align_dat_t& in_align,
             const unsigned align_no,
             const unsigned base_no,
             const RATE_GTOR_MODEL::index_t rgm,
             const bool is_codon_repeat){

  const SITE_MODEL::index_t sm(RATE_GTOR_MODEL::convert_to_site_model(rgm));
  const int offset(RATE_GTOR_MODEL::base_size_repeat_offset(rgm));

  const unsigned n_seqs(in_align.seq.dim1());

  site_code sc(n_seqs);

  for(unsigned i(0);i<n_seqs;++i){
    sc.set_taxid(i,SITE_MODEL::encode_nuc(sm,in_align.seq[i]+base_no+offset));
  }

  store_site_info(sd,in_align,align_no,base_no,sc,is_codon_repeat);
}



static
void
process_4fold_site(site_data& sd,
                   const align_dat_t& in_align,
                   const unsigned align_no,
                   const unsigned base_no){

  const unsigned n_seqs(in_align.seq.dim1());

  bool is_4fold(true);
  for(unsigned t(0);t<n_seqs;++t){
    is_4fold = is_4fold &&
      is_4fold_codon(in_align.seq[t][base_no],in_align.seq[t][base_no+1]);
  }

  if(is_4fold) {
    process_site(sd,in_align,align_no,base_no+2,RATE_GTOR_MODEL::NUC,false);
  }
}



void
nuc_seq_to_site_data(site_data& sd,
                     const nuc_seq_data& in_seq,
                     const RATE_GTOR_MODEL::index_t model,
                     const bool is_codon_border,
                     const bool is_nuc_4fold_filter,
                     const double gc_max,
                     const double gc_min,
                     const bool is_no_adjacent_nuc_diff,
                     const bool is_single_nuc_diff){

  const bool is_filter_gc(gc_max < 1. || gc_min >= 0.);

  std::vector<bool> include_seq;
  if(is_filter_gc){
    make_gc_include_map(include_seq,in_seq,gc_max,gc_min);
  }

  const unsigned n_alignments(in_seq.dat.size());
  const unsigned n_seqs(in_seq.n_seqs());

  const unsigned repeat_num(RATE_GTOR_MODEL::base_size_conditioned(model));

  bool is_codon_repeat(false);
  if      (repeat_num == 1)                { is_codon_repeat=false; }
  else if (repeat_num == CODON::BASE_SIZE) { is_codon_repeat=true; }
  else { die("unknown site model conditioned base size"); }

  const bool is_4fold(model == RATE_GTOR_MODEL::NUC && is_nuc_4fold_filter);

  sd.clear();

  for(unsigned i(0);i<n_alignments;++i){

    const align_dat_t& ali(in_seq.dat[i]);

    sd.group_label.push_back(ali.label);

    // don't accept any alignments with missing sequences:
    const unsigned alignment_n_seqs(ali.seq.dim1());
    if(n_seqs != alignment_n_seqs) continue;

    if(is_filter_gc && (! include_seq[i])) continue;

    const unsigned n_bases(ali.info.size());

    bool is_model_site_continuous(false);

    for(unsigned base_no(0);base_no<n_bases;){

      if(is_codon_repeat || is_4fold){
        // check for an actual codon:
        const bool is_codon(((base_no + (CODON::BASE_SIZE-1))<n_bases) &&
                            ali.info[base_no].codon_pos == 0 &&
                            ali.info[base_no+1].codon_pos == 1 &&
                            ali.info[base_no+2].codon_pos == 2 &&
                            ali.info[base_no+1].is_continuous &&
                            ali.info[base_no+2].is_continuous);

        if(! is_codon) {
          is_model_site_continuous = false;
          base_no++;
          continue;
        }

        const bool is_multicodon_read(is_codon_border ||
                                      RATE_GTOR_MODEL::base_size(model) > CODON::BASE_SIZE);

        if(is_multicodon_read){
          const bool is_multicodon((base_no>1) &&
                                   ((base_no+CODON::BASE_SIZE)<n_bases) &&
                                   ali.info[base_no-2].codon_pos == 1 &&
                                   ali.info[base_no-1].codon_pos == 2 &&
                                   ali.info[base_no+CODON::BASE_SIZE].codon_pos == 0 &&
                                   ali.info[base_no+CODON::BASE_SIZE+1].codon_pos == 1 &&
                                   ali.info[base_no-1].is_continuous &&
                                   ali.info[base_no].is_continuous &&
                                   ali.info[base_no+CODON::BASE_SIZE].is_continuous &&
                                   ali.info[base_no+CODON::BASE_SIZE+1].is_continuous);

          if(!is_multicodon) {
            is_model_site_continuous = false;
            base_no++;
            continue;
          }

          if(! ali.info[base_no-2].is_continuous) is_model_site_continuous = false;
        } else {
          if(! ali.info[base_no].is_continuous) is_model_site_continuous = false;
        }

      } else {
        const bool is_multinuc_read(RATE_GTOR_MODEL::base_size(model) > 1);

        if(is_multinuc_read){
          const bool is_multinuc((base_no>1) &&
                                 ali.info[base_no-1].is_continuous &&
                                 ali.info[base_no].is_continuous);

          if(!is_multinuc) {
            is_model_site_continuous = false;
            base_no++;
            continue;
          }

          if(! ali.info[base_no-2].is_continuous) is_model_site_continuous = false;
        } else {
          if(! ali.info[base_no].is_continuous) is_model_site_continuous = false;
        }
      }

      if(is_no_adjacent_nuc_diff || is_single_nuc_diff){
        bool is_bail(false);

        // only consider sites where one/zero nucleotides are different between all sites
        const SITE_MODEL::index_t sm(RATE_GTOR_MODEL::convert_to_site_model(model));
        const int offset(RATE_GTOR_MODEL::base_size_repeat_offset(model));
        const int size(SITE_MODEL::base_size(sm));

        const int start(static_cast<int>(base_no)+offset);
        if(start>=0 && (start+size)<static_cast<int>(n_bases)){
          unsigned n_diff(0);
          bool is_last_diff(false);
          for(int bb(start);bb<(start+size);++bb){
            bool is_diff(false);
            const NUC::index_t n0(ali.seq[0][bb]);
            for(unsigned j(1);j<alignment_n_seqs;++j){
              if(n0 != ali.seq[j][bb]){
                n_diff++;
                is_diff=true;
                break;
              }
            }
            if(is_no_adjacent_nuc_diff && is_diff && is_last_diff) {
              is_bail=true;
              break;
            }
            is_last_diff=is_diff;
          }
          if(is_single_nuc_diff) is_bail=(n_diff>1);
        } else { is_bail=true; }

        if(is_bail){
          is_model_site_continuous = false;
          base_no++;
          continue;
        }
      }

      if(is_4fold){
        process_4fold_site(sd,ali,i,base_no);
        base_no += CODON::BASE_SIZE;
      }else {
        process_site(sd,ali,i,base_no,model,is_codon_repeat);
        base_no += repeat_num;
      }

      is_model_site_continuous = true;
    }
  }

  sd.sm = RATE_GTOR_MODEL::convert_to_site_model(model);
  sd.data_class_labels = in_seq.data_class_labels;
  sd.taxid = in_seq.taxid;
}

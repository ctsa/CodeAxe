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

// $Id: process_fsa.cc 1203 2008-05-05 19:17:32Z ctsa $

/// \file

#include "process_fsa.h"
#include "simple_util.h"
#include "util/bio/fasta_seq.h"
#include "util/bio/seq_util.h"
#include "util/general/die.h"
#include "util/general/log.h"
#include "util/general/io_util.h"

#include <sys/types.h>
#include <dirent.h>

#include <cassert>
#include <cctype>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using namespace std;


const char DEFAULT_CLASS_LABEL[] = "none";

const unsigned CODON_BASE_SIZE(3);

typedef std::vector<simple_init_array<char> > tmp_seq_type;


// convert input sequence to align_dat_t format, no
// questions asked:
//
static
void
no_filter(align_dat_t& filter_align,
          const align_dat_t::info_type& info,
          const tmp_seq_type& seq){

  const unsigned n_seqs(seq.size());
  if(n_seqs==0) return;

  const unsigned n_bases(seq[0].size());

  filter_align.info.resize(n_bases);
  filter_align.seq.init(n_seqs,n_bases);

  for(unsigned i(0);i<n_bases;++i){
    filter_align.info[i] = info[i];
    for(unsigned s(0);s<n_seqs;++s){
      filter_align.seq[s][i] = NUC::remap(seq[s][i]);
    }
  }
}



// check for 'check-val' in any seq at codon site::
static
bool
check_codon_site(const tmp_seq_type& seq,
                 const unsigned codonno,
                 const char check_val){

  const unsigned baseno(codonno*CODON_BASE_SIZE);
  const unsigned n_seqs(seq.size());

  for(unsigned t(0);t<n_seqs;++t){
    const char aa_site(codon_trans_known(seq[t].ptr()+baseno));
    if (aa_site == check_val) { return true; }
  }
  return false;
}


// check for stop in current position:
static
bool
check_stop_codon(const tmp_seq_type& seq,
                 const unsigned codonno){

  return check_codon_site(seq,codonno,'*');
}


// check for gap/ambiguity/stop in current position::
static
bool
check_gap_codon(const tmp_seq_type& seq,
                const unsigned codonno){

  return check_codon_site(seq,codonno,'-');
}


static
bool
check_ambi_codon(const tmp_seq_type& seq,
                 const unsigned codonno){

  return check_codon_site(seq,codonno,'X');
}




// assume this is aligned coding sequence, skip gaps and codons with
// ambiguities, previous and next codons must be clean as well,
// gapless filter currently applies to the entire column only, which
// only works as long as we allow no missing sequences/codons
//
static
void
codon_filter(align_dat_t& filter_align,
             const align_dat_t::info_type& info,
             const std::string& label, //only required for error report
             const tmp_seq_type& seq,
             const seq_filter_options& fopt){

  static const unsigned default_seqwindow(10);
  static const unsigned min_gapless(20);

  const unsigned n_seqs(seq.size());
  if(n_seqs==0) return;

  const unsigned n_bases(seq[0].size());

  const unsigned seqwindow((fopt.is_use_window_filter?default_seqwindow:0));

  if( n_bases%CODON_BASE_SIZE ) {
    std::ostringstream oss;
    oss << "process_alignment codon_filter() fsa sequences not divisible to codons, fsa label: " << label;

    die(oss.str().c_str());
  }

  const unsigned n_codons(n_bases/CODON_BASE_SIZE);

  std::vector<bool> pass_list(n_codons,false);

  // pre-calc seq aa's
  simple_matrix<char> aa(n_seqs,n_codons);
  for(unsigned t(0);t<n_seqs;++t){
    for(unsigned codonno=0;codonno<n_codons;codonno++){
      const unsigned baseno(codonno*CODON_BASE_SIZE);
      aa[t][codonno] = codon_trans(seq[t].ptr()+baseno);
    }
  }

  // go thru entire sequence, bounded by window-ends...
  //
  for(unsigned codonno(seqwindow); (codonno+seqwindow) < n_codons; ++codonno){

    if(check_stop_codon(seq,codonno)) continue;
    if(fopt.is_use_gap_filter && check_gap_codon(seq,codonno)) continue;
    if(fopt.is_use_ambi_filter && check_ambi_codon(seq,codonno)) continue;

    // run windowed quality filter:
    if(fopt.is_use_window_filter){

      unsigned n_gapless(0);
      unsigned n_match(0);

      for(int i=-static_cast<int>(seqwindow);i<(static_cast<int>(seqwindow)+1);++i){
        if (i==0) continue;

        const unsigned window_codonno(codonno+i);

        bool is_gap(false);
        for(unsigned t(0);t<n_seqs;++t){
          if(aa[t][window_codonno] == '-'){
            is_gap=true;
            break;
          }
        }

        if(! is_gap) n_gapless++;

        bool is_match(true);
        char maa(0);
        for(unsigned t(0);t<n_seqs;++t){
          if(t==0){
            maa = aa[t][window_codonno];
          } else {
            if( maa != aa[t][window_codonno] ){
              is_match=false;
              break;
            }
          }
        }

        if(is_match) n_match++;
      }

      // filter from window results
      if(n_gapless < min_gapless) continue;
      if(n_match < fopt.min_match) continue;
    }

    // made it! record codon:
    pass_list[codonno] = true;
  }

  unsigned n_pass_bases(0);
  for(unsigned i(0);i<n_codons;++i) if(pass_list[i]) n_pass_bases+=3;

  filter_align.info.resize(n_pass_bases);
  filter_align.seq.init(n_seqs,n_pass_bases);

  bool is_last_accepted_codonno_set(false);
  unsigned last_accepted_codonno(0);
  for(unsigned i(0),pi(0);i<n_codons;++i){
    const bool is_continuous(is_last_accepted_codonno_set && (i == last_accepted_codonno+1));

    if(pass_list[i]){
      for(unsigned j(0);j<CODON_BASE_SIZE;++j){
        const unsigned base(i*CODON_BASE_SIZE+j);
        const unsigned pbase(pi*CODON_BASE_SIZE+j);
        filter_align.info[pbase] = info[base];
        filter_align.info[pbase].is_continuous = (j!=0 || (filter_align.info[pbase].is_continuous && is_continuous));
        for(unsigned s(0);s<n_seqs;++s){
          filter_align.seq[s][pbase] = NUC::remap(seq[s][base]);
        }
      }
      pi++;

      is_last_accepted_codonno_set = true;
      last_accepted_codonno = i;
    }
  }
}



static
void
harvest_fsa_file(align_dat_t::info_type& info,
                 tmp_seq_type& seq,
                 name_id_lup& taxid_lup,
                 const char* fsa_file){

  // just assume it's the file for now
  ifstream infp(fsa_file);
  check_nonempty_istream(infp,fsa_file);

  info.clear();
  seq.clear();

#if 0
  // faster code to replace FastaSeq functions
  //
  int c=infp.peek();
  if(c!='>') die("can't parse fasta");

  // read seq label from the 3rd column:
  string seq_label;
  infp >> seq_label >> seq_label >> seq_label;

  // kill ws:
  do{
    c=infp.get();
  } while(c==' ' || c=='\t' || c== '\n' || ! infp);

  infp.unget(c);

  tseq.resize(0);

  while(true){
    c=infp.get();
    if(! infp) break;
    if(c==' ' || c=='\t' || c== '\n') continue;
    if(c=='>') {
      infp.unget(c);
      break;
    }
    tseq.push_back(c);
  }

  tseq.push_back
#endif

  FastaSeq f;

  unsigned seqcount(0);

  while( infp >> f ) {
    seqcount++;

    istringstream is(f.header);

    // throw away the '>' char on header
    assert(is.peek() == '>');
    is.get();

    // first ws separated word in header is used as org name
    string taxstr;
    is >> taxstr;

    const unsigned taxid(taxid_lup.assignid(taxstr));
    if(taxid>=seq.size()){ seq.resize(taxid+1); }
    seq[taxid].init(f.seq.size());
    std::copy(f.seq.begin(),f.seq.end(),seq[taxid].ptr());
  }

  const unsigned n_seqs(seq.size());

  // check that all seq are the same size:
  for(unsigned t(0),n_bases(0);t<n_seqs;++t){
    const unsigned b(seq[t].size());
    if(n_bases==0) n_bases = b;
    if(b==0){
      // sequence absent: non-error skip case
      seq.clear();
      return;
    } else if(n_bases != b){
      die("process_fsa:: invalid alignment format");
    }
  }

  if(n_seqs){
    const unsigned n_bases(seq[0].size());

    // add default info
    info.resize(n_bases);
    for(unsigned n(0);n<n_bases;++n){
      info[n].class_no = DEFAULT_CLASS_ID;
      info[n].codon_pos = static_cast<unsigned char>(n%CODON_BASE_SIZE);
      info[n].is_continuous = (n!=0);
    }

    // upper-caseify everything:
    for(unsigned t(0);t<n_seqs;++t){
      char* v(seq[t].ptr());
      for(unsigned n(0);n<n_bases;++n){ v[n] = toupper(v[n]); }
    }
  }
}




// assumes in_align is sorted on base_no!!
//
// parsing is lazy sstream stuff, should not scale well.
//
// assumes input data classes are numbered by codon (zero-indexed)
//
//format:
// first line: [codon_no|*] class_no
// other line: codon_no class_no
//
// * sets the default class_no for all valid codons in the sequence
// else, default class_no=0
//
static
void
harvest_class_file(align_dat_t::info_type& info,
                   const char* class_file){

  // class file is allowed to be empty or not exist to indicate default class:
  ifstream class_is(class_file);
  if( ! class_is ) return;

  string buf;
  bool is_read_first_noncomment(false);

  const unsigned asize(info.size());

  unsigned file_default_class_id(DEFAULT_CLASS_ID);

  while(getline(class_is,buf)){
    if( buf[0] == '#') continue;

    istringstream sbuf(buf);

    unsigned codon_no;
    unsigned class_no;

    if( ! is_read_first_noncomment ){
      is_read_first_noncomment=true;
      if(buf[0] == '*'){
        sbuf.get();
        sbuf >> class_no;
        if(!sbuf) die((string("irregular class file format (1): ")+class_file).c_str());

        if(class_no == file_default_class_id) continue;

        for(unsigned i(0);i<asize;++i){
          info[i].class_no = static_cast<align_col_t::class_type>(class_no);
        }

        // after setting the default, individual codon data classes can be
        // specified, this lets downstream know that a new default has
        // been set for this alignment
        file_default_class_id = class_no;
        continue;
      }
    }

    sbuf >> codon_no >> class_no;

    if(!sbuf) die((string("irregular class file format (2): ")+class_file).c_str());

    // skip default class value:
    if(class_no == file_default_class_id) continue;

    const unsigned base_no(codon_no*CODON_BASE_SIZE);
    for(unsigned j(0);j<CODON_BASE_SIZE;++j) {
      info[base_no+j].class_no = static_cast<align_col_t::class_type>(class_no);
    }
  }
}



/// class_label format:
/// # comment anyline
/// 1 label1
/// 2 label2
/// ...
///
static
void
harvest_class_label(map<unsigned,string>& id_class_label_map,
                    const char* class_label_file){

  string buf;
  ifstream class_is(class_label_file);
  check_nonempty_istream(class_is,class_label_file);

  while(getline(class_is,buf)){
    if( buf[0] == '#' ) continue;
    istringstream sbuf(buf);

    unsigned class_id;
    string class_label;
    sbuf >> class_id >> class_label;

    if(class_label.size() == 0) pass_away("empty label in class_label file");

    if(id_class_label_map.find(class_id) != id_class_label_map.end()){
      if(id_class_label_map[class_id] != class_label &&
         ! (class_id == DEFAULT_CLASS_ID && id_class_label_map[class_id] == DEFAULT_CLASS_LABEL)){
        log_os << "ERROR:: harvest_class_label: class_od: " << class_id
               << " is already assigned to label: " << id_class_label_map[class_id] << "\n";
        exit(EXIT_FAILURE);
      }
    }
    id_class_label_map[class_id] = class_label;
  }
}




static
void
harvest_file_set(nuc_seq_data& in_seq,
                 align_dat_t::info_type& tmp_info,
                 tmp_seq_type& tmp_seq,
                 const char* full_fsa_file,
                 const bool is_data_class_mode,
                 const char* full_class_file,
                 const seq_filter_options& fopt,
                 const std::string& label){

  harvest_fsa_file(tmp_info,tmp_seq,in_seq.taxid,full_fsa_file);
  if(tmp_seq.empty()) return;

  if(is_data_class_mode){
    harvest_class_file(tmp_info,full_class_file);
  }

  in_seq.dat.push_back(align_dat_t());
  if(fopt.is_skip_filters){
    no_filter(in_seq.dat.back(),tmp_info,tmp_seq);
  } else {
    codon_filter(in_seq.dat.back(),tmp_info,label,tmp_seq,fopt);
  }
  in_seq.dat.back().label = label;
}




const char dir_delim('/');
const char* fsa_ext(".fsa");
const char* class_ext(".class");




void
process_seq_data(nuc_seq_data& in_seq,
                 const seq_filter_options& fopt,
                 const char* seq_file_node,
                 const char* class_file_node,
                 const char* class_labels_file){

  if(! (seq_file_node && strlen(seq_file_node))){
    die("invalid input sequence file/directory");
  }

  // initialize in_seq:
  in_seq.clear();

  std::map<unsigned,std::string> id_class_label_map;
  id_class_label_map[DEFAULT_CLASS_ID]=DEFAULT_CLASS_LABEL;

  bool is_data_class_mode(false);
  if(class_file_node && strlen(class_file_node)){
    is_data_class_mode=true;
    if(class_labels_file && strlen(class_labels_file)){
      harvest_class_label(id_class_label_map,class_labels_file);
    }
  }

  DIR *dir(opendir(seq_file_node));
  DIR *class_dir(0);
  if(is_data_class_mode) class_dir=opendir(class_file_node);

  align_dat_t::info_type tmp_info;
  tmp_seq_type tmp_seq;

  if(dir){
    log_os << "reading sequence directory: " << seq_file_node << "\n";
    if(is_data_class_mode){
      if(class_dir==0) die("expected data class directory");
      log_os << "      data class directory: " << class_file_node << "\n";
      closedir(class_dir);
    }

    string full_fsa_file;
    string full_class_file;
    unsigned attempted_align_count(0);

    dirent* ent;
    while (0!=(ent=readdir(dir))){
      if( strstr(ent->d_name,fsa_ext) == NULL ) continue;
      const char* fsa_file(ent->d_name);

      full_fsa_file=string(seq_file_node)+dir_delim+fsa_file;

      // strip off filename suffix and store label:
      string fsa_label(fsa_file);
      size_t a(fsa_label.rfind('.'));
      if(a != string::npos) fsa_label.erase(a);

      if(is_data_class_mode){
        full_class_file=string(class_file_node)+dir_delim+fsa_label+class_ext;
      }

      harvest_file_set(in_seq,tmp_info,tmp_seq,full_fsa_file.c_str(),is_data_class_mode,
                       full_class_file.c_str(),fopt,fsa_label);

      attempted_align_count++;
      const unsigned readable_align_count(in_seq.dat.size());

      if(attempted_align_count%20==0) log_os << ".";
      if(attempted_align_count%1000==0) log_os << "\n";
      if(fopt.max_align_count>0 && readable_align_count>=fopt.max_align_count) break;
    }
    closedir(dir);
    log_os << "\n";

  } else {
    log_os << "reading sequence file: " << seq_file_node << "\n";

    string fsa_label(seq_file_node);
    size_t a(fsa_label.rfind('.'));
    if(a != string::npos) fsa_label.erase(a);

    if(is_data_class_mode){
      if(class_dir) die("expected data class file");
      log_os << "        data class file: " << class_file_node << "\n";
    }

    harvest_file_set(in_seq,tmp_info,tmp_seq,seq_file_node,is_data_class_mode,
                     class_file_node,fopt,fsa_label);
  }

  { // make compressed id-label map for storage of data:
    std::map<unsigned,std::string>::const_iterator i,i_end=id_class_label_map.end();
    for(i=id_class_label_map.begin();i!=i_end;++i) {
      in_seq.data_class_labels.assignid(i->second);
    }
  }

  // convert class_ids to compressed values:
  {
    std::map<unsigned,unsigned> conv_map;
    std::map<unsigned,std::string>::const_iterator i,i_end=id_class_label_map.end();
    for(i=id_class_label_map.begin();i!=i_end;++i) {
      conv_map[i->first] = in_seq.data_class_labels.getid(i->second);
    }

    const unsigned n_alignments(in_seq.dat.size());
    for(unsigned j(0);j<n_alignments;++j){
      align_dat_t& ad(in_seq.dat[j]);
      const unsigned n_cols(ad.info.size());
      for(unsigned k(0);k<n_cols;++k){
        ad.info[k].class_no = conv_map[ad.info[k].class_no];
      }
    }
  }
}

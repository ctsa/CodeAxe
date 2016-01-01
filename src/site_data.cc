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

// $Id: site_data.cc 1193 2008-03-29 03:14:54Z ctsa $

/// \file

#include "site_data.h"
#include "util/bio/bioseq_util.h"
#include "util/general/die.h"
#include "util/general/io_util.h"
#include "util/general/log.h"
#include "util/general/metatags.h"

#include <cassert>
#include <cstdlib>

#include <iostream>
#include <sstream>
#include <vector>
#include <utility>

using namespace std;


const unsigned DATA_FILE_VERSION(3);



void
site_data::
peek(ostream& os) const {
  static const unsigned peek_len(10);

  site_count_map::const_iterator i,i_end=site_count.end();

  const unsigned len(site_count.size());
  os << "seq_t:\n";
  os << "len: " << len << "\n";

  unsigned j(0);
  for(i=site_count.begin();i!=i_end;++i){
    if(j++>=peek_len) break;
    os << i->first << "\n";
  }
}



static
void
finish_data_class_line(const vector<unsigned>& line_class_count,
                       const unsigned n_line_classes,
                       ostream& os){
  for(unsigned c(0),sum(0);c<n_line_classes;++c){
    // check for and skip printing of final "1" if it's the only count
    // in this line (@ 1group/line):
    if((c+1)==n_line_classes && line_class_count[c]==1 && sum == 0) continue;
    os << " " << line_class_count[c];
    sum += line_class_count[c];
  }
  os << "\n";
}



static const char* org_signal("org:");
static const char* dc_signal("data_class:");
static const char* tmp_cat_signal("cat:");
static const char* glabel_signal("group_labels:");
static const char* sitedata_signal("sites:");
static const char site_signal('S');




void
site_data::
store_state(ostream& os) const {

  os << "site_data_file_version " << DATA_FILE_VERSION << "\n";
  _ai.store_state(os);
  os << "site_model " << sm << "\n";

  const unsigned n_org(taxid.size());

  os << org_signal;
  for(unsigned t(0);t<n_org;++t){
    os << " " << taxid.getstr(t);
  }
  os << "\n";

  const unsigned n_classes(data_class_size());

  os << dc_signal;
  for(unsigned i(0);i<n_classes;++i){
    os << " " << data_class_labels.getstr(i);
  }
  os << "\n";

  os << glabel_signal << "\n";
  const unsigned gsi(group_label.size());
  for(unsigned i(0);i<gsi;++i){
    os << group_label[i] << "\n";
  }

  os << sitedata_signal << "\n";

  vector<unsigned> line_class_count(n_classes);

  site_count_map::const_iterator i,i_end=site_count.end();

  for(i=site_count.begin();i!=i_end;++i){
    os << site_signal << i->first;

    bool is_new_group(true);
    unsigned last_group(0);
    unsigned n_line_classes(0);

    site_count_data_type::const_iterator j,j_end=i->second.end();

    for(j=i->second.begin();j!=j_end;++j) {
      const site_data_code::group_id_type&  group(j->first.group_id);
      const site_data_code::data_class_id_type& class_id(j->first.data_class_id);
      const unsigned& count(j->second);

      if(group != last_group) is_new_group=true;

      if(is_new_group) {
        finish_data_class_line(line_class_count,n_line_classes,os);
        os << group-last_group;
        n_line_classes=0;
        is_new_group=false;
      }

      const unsigned last_nlc(n_line_classes);
      n_line_classes=class_id+1;

      for(unsigned c(last_nlc);c<n_line_classes;++c) {
        if   (c==class_id) line_class_count[c] = count;
        else               line_class_count[c] = 0;
      }
      last_group = group;
    }
    finish_data_class_line(line_class_count,n_line_classes,os);
  }
}



static
void
get_site_code(const unsigned nseq,
              site_code& sc,
              istream& is){

  int c;
  unsigned sci(0);
  site_code::index_type it;

  while( (c=is.get()) != '\n' && is && sci < nseq){
    if(c==' ' || c=='\t') continue;
    is.unget();
    is >> it;
    sc.set_taxid(sci,it);
    sci++;
  }
}



static
void
load_die(istream& is,
         const char* const sid,
         const string& tmpstr)  NORETURN_TAG;

static
void
load_die(istream&,
         const char* const sid,
         const string& tmpstr){

  log_os << "ERROR:: -- Invalid site data file format. --\n"
         << "ERROR:: expected_section: " << sid << "\n"
         << "ERROR:: section: " << tmpstr << "\n"
         << "\n";
  //  log_os << is.rdbuf();
  abort();
}




static
void
load_advance(istream& is,
             const char* const second = 0){

  string tmpstr;
  is >> tmpstr;
  if(second && tmpstr != second) load_die(is,second,tmpstr);
}




void
site_data::
load_state(istream& is,
           const bool is_read_group_info) {

  clear();

  {  // check file version
    unsigned tmpu;
    load_advance(is,"site_data_file_version");
    is >> tmpu;
    if(tmpu != DATA_FILE_VERSION){
      log_os << "ERROR:: bad site data file version number. expected: " << DATA_FILE_VERSION << "\n";
      exit(EXIT_FAILURE);
    }
  }

  // throw away old audit_info
  _ai.skip_state(is);

  // get site model
  {
    string dummy;
    is >> dummy;
    sm = static_cast<SITE_MODEL::index_t>(read_int(is));
  }
  clear_line(is);

  {
    string buf;

    while(is){
      if(is.peek()==site_signal) break;

      if( ! getline(is,buf) ) break;

      istringstream sbuf(buf);
      string tmpstr;

      sbuf >> tmpstr;
      if       (tmpstr==org_signal){
        while(sbuf >> tmpstr){ taxid.assignid(tmpstr); }
      } else if( tmpstr==dc_signal || tmpstr==tmp_cat_signal){ // backcompatibility
        while(sbuf >> tmpstr) { data_class_labels.assignid(tmpstr); }
      } else if(tmpstr==glabel_signal){
        break;
      } else {
        die((string("unexpected site_data format at line:\n")+buf).c_str());
      }
    }

    while(is){
      if( ! getline(is,buf) ) break;
      if(buf==sitedata_signal) break;

      group_label.push_back(buf);
    }
  }

  const unsigned nseq(taxid.size());
  if(nseq==0) die("no sequence labels in site data file");

  site_code sc(nseq);

  int c;

  unsigned group_no(0);
  unsigned tmp_uint;
  unsigned data_class_id;
  bool is_zero_line;

  site_count_data_type* site_map(0);

  while(is){
    c=is.get();
    if(c==' ' || c=='\t' || c== '\n' || ! is) continue;
    if(c==site_signal){

      group_no=0;
      get_site_code(nseq,sc,is);

      site_map = &(site_count[sc]);

    } else {
      is.unget();
      is >> tmp_uint; // group_delta
      if(is_read_group_info) group_no += tmp_uint;
      data_class_id=0;
      is_zero_line=true;
      while( (c=is.get()) != '\n' && is){
        if(c==' ' || c=='\t') continue;
        is.unget();
        is >> tmp_uint; // count

        if(tmp_uint != 0) {
          site_map->operator[](site_data_code(group_no,data_class_id)) += tmp_uint;
          is_zero_line=false;
        }
        data_class_id++;
      }
      if(is_zero_line) { site_map->operator[](site_data_code(group_no,data_class_id)) += 1; }
    }
  }
}

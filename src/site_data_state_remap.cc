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

// $Id: site_data_state_remap.cc 1163 2008-03-21 00:13:32Z ctsa $

/// \file

#include "site_data.h"
#include "site_data_state_remap.h"
#include "util/bio/bioseq_util.h"
#include "util/math/prob_util.h"

using namespace std;

template <typename StateCompressor>
void
site_data_state_compressor(site_data& sd,
                           const StateCompressor& scom){

  site_count_map::const_iterator i,i_end;

  site_count_map site_count_out;

  i=sd.site_count.begin();
  i_end=sd.site_count.end();

  const unsigned n_seqs(sd.n_taxa());
  site_code sc(n_seqs);

  for(;i!=i_end;++i){
    for(unsigned t(0);t<n_seqs;++t){
      sc.set_taxid(t,scom(i->first.get_taxid(t)));
    }

    site_count_data_type::const_iterator j,j_end;

    j=i->second.begin();
    j_end=i->second.end();

    for(;j!=j_end;++j){
      site_count_out[sc][j->first] += j->second;
    }
  }

  sd.site_count = site_count_out;
}


struct nsc5_to_nsc4pre {
  unsigned operator()(unsigned x) const {
    NSCODON::index_t c;
    NUC::index_t n_pre,n_post;
    NSC5::decode(n_pre,n_post,c,static_cast<NSC5::index_t>(x));
    return NSC4::encode(n_pre,c);
  }
};

struct nsc5_to_nsc4post {
  unsigned operator()(unsigned x) const {
    NSCODON::index_t c;
    NUC::index_t n_pre,n_post;
    NSC5::decode(n_pre,n_post,c,static_cast<NSC5::index_t>(x));
    return NSC4::encode(n_post,c);
  }
};

struct nsc5_to_nscodon {
  unsigned operator()(unsigned x) const {
    return NSC5::decode_nscodon(static_cast<NSC5::index_t>(x));
  }
};

struct nsc4_to_nscodon {
  unsigned operator()(unsigned x) const {
    return NSC4::decode_nscodon(static_cast<NSC4::index_t>(x));
  }
};


void
site_data_nsc5_to_nsc4pre(site_data& sd){
  site_data_state_compressor(sd,nsc5_to_nsc4pre());
}


void
site_data_nsc5_to_nsc4post(site_data& sd){
  site_data_state_compressor(sd,nsc5_to_nsc4post());
}


void
site_data_nsc5_to_nscodon(site_data& sd){
  site_data_state_compressor(sd,nsc5_to_nscodon());
}


void
site_data_nsc4_to_nscodon(site_data& sd){
  site_data_state_compressor(sd,nsc4_to_nscodon());
}

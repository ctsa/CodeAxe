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

// $Id: site_data_stats.cc 1170 2008-03-24 20:55:56Z ctsa $

/// \file

#include "site_data.h"
#include "site_data_stats.h"
#include "util/bio/bioseq_util.h"
#include "util/math/prob_util.h"

#include <cstring>

using namespace std;


// get expected symbol distribution from all site data
void
site_data_2_state_pdf(const site_data& sd,
                      const unsigned n_states,
                      prob_t* distro){

  const unsigned n_seqs(sd.n_taxa());

  memset(distro,0,n_states*sizeof(prob_t));

  site_count_map::const_iterator i,i_end=sd.site_count.end();

  for(i=sd.site_count.begin();i!=i_end;++i){
    site_count_data_type::const_iterator j,j_end=i->second.end();

    for(j=i->second.begin();j!=j_end;++j){
      for(unsigned t(0);t<n_seqs;++t){
        const unsigned state(i->first.get_taxid(t));
        if(state>=n_states) continue;
        distro[state] += static_cast<prob_t>(j->second);
      }
    }
  }

  pdistro_norm(distro,distro+n_states);
}


template <typename T>
void
site_data_2_taxa_state_count_internal(const site_data& sd,
                                      const unsigned n_states,
                                      T* count,
                                      const bool is_accumulate = false,
                                      const bool is_single_data_class = false,
                                      const unsigned target_data_class_id = 0){

  const unsigned n_seqs(sd.n_taxa());

  if(! is_accumulate) std::fill(count,count+n_states*n_seqs,static_cast<T>(0));

  site_count_map::const_iterator i,i_end=sd.site_count.end();

  for(i=sd.site_count.begin();i!=i_end;++i){
    site_count_data_type::const_iterator j,j_end=i->second.end();
    for(j=i->second.begin();j!=j_end;++j){
      if(is_single_data_class){
        const unsigned data_class_id(static_cast<unsigned>(j->first.data_class_id));
        if(data_class_id != target_data_class_id) continue;
      }
      for(unsigned t(0);t<n_seqs;++t){
        const unsigned state(i->first.get_taxid(t));
        if(state>=n_states) continue;  // ambiguous state

        count[t*n_states+state] += static_cast<T>(j->second);
      }
    }
  }
}


// get state counts for each org from site data
void
site_data_2_taxa_state_count(const site_data& sd,
                             const unsigned n_states,
                             unsigned* count,
                             const bool is_accumulate,
                             const bool is_single_data_class,
                             const unsigned target_data_class_id){

  site_data_2_taxa_state_count_internal(sd,n_states,count,is_accumulate,is_single_data_class,target_data_class_id);
}


// get state distribution for each org from site data
void
site_data_2_taxa_state_pdf(const site_data& sd,
                           const unsigned n_states,
                           prob_t* distro){

  const unsigned n_seqs(sd.n_taxa());

  site_data_2_taxa_state_count_internal(sd,n_states,distro);

  for(unsigned t(0);t<n_seqs;++t) {
    pdistro_norm(distro+t*n_states,distro+(t+1)*n_states);
  }
}



template <typename RandomAccessIterator>
void
site_data_2_data_class_state_count_internal(const site_data& sd,
                                            const unsigned n_data_classes,
                                            const unsigned n_states,
                                            RandomAccessIterator data_class_state_count){

  const unsigned n_seqs(sd.n_taxa());

  std::fill(data_class_state_count,data_class_state_count+n_data_classes*n_states,0);

  site_count_map::const_iterator i,i_end=sd.site_count.end();

  for(i=sd.site_count.begin();i!=i_end;++i){
    site_count_data_type::const_iterator j,j_end=i->second.end();

    for(j=i->second.begin();j!=j_end;++j){
      const unsigned data_class_id(static_cast<unsigned>(j->first.data_class_id));
      if(data_class_id>=n_data_classes) continue;
      for(unsigned t(0);t<n_seqs;++t){
        const unsigned state(i->first.get_taxid(t));
        if(state>=n_states) continue;
        *(data_class_state_count+data_class_id*n_states+state) += j->second;
      }
    }
  }
}



void
site_data_2_data_class_state_count(const site_data& sd,
                                   const unsigned n_data_classes,
                                   const unsigned n_states,
                                   std::vector<unsigned>& data_class_state_count){

  data_class_state_count.resize(n_data_classes*n_states);
  site_data_2_data_class_state_count_internal(sd,n_data_classes,n_states,
                                              data_class_state_count.begin());
}



void
site_data_2_data_class_state_pdf(const site_data& sd,
                                 const unsigned n_data_classes,
                                 const unsigned n_states,
                                 prob_t* distro) {

  site_data_2_data_class_state_count_internal(sd,n_data_classes,n_states,distro);

  for(unsigned c(0);c<n_data_classes;++c){
    pdistro_norm_safe(distro+c*n_states,distro+(c+1)*n_states);
  }
}



void
site_data_2_data_class_count(const site_data& sd,
                             const unsigned n_data_classes,
                             std::vector<unsigned>& data_class_count){

  site_count_map::const_iterator i,i_end=sd.site_count.end();

  for(i=sd.site_count.begin();i!=i_end;++i){
    site_count_data_type::const_iterator j,j_end=i->second.end();

    for(j=i->second.begin();j!=j_end;++j){
      const unsigned data_class_id(static_cast<unsigned>(j->first.data_class_id));
      if(data_class_id>=n_data_classes) continue;
      data_class_count[data_class_id] += j->second;
    }
  }
}



void
site_data_counts(const site_data& sd,
                 unsigned& site_count,
                 unsigned& conserved_site_count,
                 unsigned& unique_site_count){

  site_count=0;
  conserved_site_count=0;
  unique_site_count=sd.site_count.size();

  const unsigned n_seqs(sd.n_taxa());

  site_count_map::const_iterator i,i_end=sd.site_count.end();

  for(i=sd.site_count.begin();i!=i_end;++i){
    site_count_data_type::const_iterator j,j_end=i->second.end();

    // check if this is a conserved site:
    bool is_conserved(true);
    {
      unsigned last_state(0);
      for(unsigned t(0);t<n_seqs;++t){
        const unsigned state(i->first.get_taxid(t));
        if(t!=0 && state!=last_state) {
          is_conserved=false;
          break;
        }
        last_state=state;
      }
    }

    // sum over all categories:
    for(j=i->second.begin();j!=j_end;++j){
      site_count += j->second;
      if(is_conserved) conserved_site_count += j->second;
    }
  }
}


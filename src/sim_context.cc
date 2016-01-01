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

// $Id: sim_context.cc 1193 2008-03-29 03:14:54Z ctsa $

/// \file

#include "bi_tree.h"
#include "cat_manager.h"
#include "cat_info.h"
#include "rate_gtor_nscodon_base.h"
#include "rate_gtor_nuc_base.h"
#include "root_gtor.h"
#include "sim_context.h"
#include "sim_context_continuous.h"
#include "sim_context_discrete.h"
#include "sim_context_shared.h"
#include "sim_util.h"
#include "simple_util.h"
#include "time_gtor.h"
#include "util/bio/bioseq_util.h"
#include "util/general/die.h"
#include "util/general/log.h"
#include "util/math/array_util.h"
#include "util/math/prob_util.h"
#include "util/math/random_util.h"

#include <ostream>
#include <utility>
#include <vector>




const unsigned NUC3MER_SIZE(NUC::SIZE*NUC::SIZE*NUC::SIZE);


// for noncoding models
//
static
void
nuc_mut_prob_context_from_single_context(const rate_gtor_nuc_base& r,
                                         const NUC::index_t n_minus1,
                                         const NUC::index_t n,
                                         const NUC::index_t n_plus1,
                                         prob_t& subprob,
                                         prob_t* subdist,
                                         const unsigned cat,
                                         const unsigned branch_cat_set){

  // loop through possible nuc changes at n:
  //
  for(unsigned j(0);j<NUC::SIZE;++j){
    const NUC::index_t nj(static_cast<NUC::index_t>(j));
    if(nj == n) continue;

    smlfloat this_rate = r.cat_nuc_context_mut_rate(n_minus1,n,n_plus1,nj,cat,branch_cat_set,0,0);

    subdist[nj+(nj<n?0:-1)] = this_rate;
  }
  subprob = array_sum(subdist,(NUC::SIZE-1));
}



// for noncoding models
//
static
void
nuc_mut_prob_context_from_all_contexts(const rate_gtor_nuc_base& r,
                                       prob_t* nucsub,
                                       prob_t* nucsub_cdf,
                                       const unsigned cat,
                                       const unsigned branch_cat_set,
                                       const bool is_calculate_cdf){

  for(unsigned n0(0);n0<NUC::SIZE;++n0){
    for(unsigned n1(0);n1<NUC::SIZE;++n1){
      for(unsigned n2(0);n2<NUC::SIZE;++n2){
        const unsigned nuc3mer(n0+n1*NUC::SIZE+n2*NUC::SIZE*NUC::SIZE);

        prob_t subdist[(NUC::SIZE-1)];
        nuc_mut_prob_context_from_single_context(r,
                                                 static_cast<NUC::index_t>(n0),
                                                 static_cast<NUC::index_t>(n1),
                                                 static_cast<NUC::index_t>(n2),
                                                 nucsub[nuc3mer],subdist,
                                                 cat,branch_cat_set);

        // load subdist into exported data structure:
        if( is_calculate_cdf ){
          prob_t* nucsub_cdf_3mer = nucsub_cdf+nuc3mer*(NUC::SIZE-1);
          for(unsigned k(0);k<(NUC::SIZE-1);++k){
            nucsub_cdf_3mer[k] = subdist[k];
          }
        }
      }
    }
  }
}



//
static
void
nuc_mut_prob_context(const rate_gtor_nuc_base& r,
                     prob_t** nucsub,       // [n_cats][NUC3MER_SIZE]
                     prob_t** nucsub_cdf,   // [n_cats][NUC3MER_SIZE*(NUC::SIZE-1)]
                     const bool is_use_discrete_time,
                     const smlfloat discrete_time_unit_normed,
                     const unsigned n_cats,
                     const time_gtor& tg,
                     const unsigned branch_no,
                     const cat_manager& cm,
                     const bool is_calculate_cdf = true){

  const smlfloat etime(tg.cat_expected_branch_time(branch_no));

  for(unsigned c(0);c<n_cats;++c){
    const unsigned bcs(cm.get_branch_cat_set(c,branch_no));

    prob_t* nucsub_cat(nucsub[c]);
    prob_t* nucsub_cdf_cat(0);
    if( is_calculate_cdf ) nucsub_cdf_cat = nucsub_cdf[c];

    nuc_mut_prob_context_from_all_contexts(r,nucsub_cat,nucsub_cdf_cat,
                                           c,bcs,is_calculate_cdf);

    const smlfloat branch_time(tg.branch_time(branch_no,c));
    const smlfloat branch_cat_scale(branch_time <= 0. ? 0. : branch_time/etime);
    array_scale(nucsub_cat,NUC3MER_SIZE,branch_cat_scale);

    if(is_use_discrete_time){
      array_scale(nucsub_cat,NUC3MER_SIZE,discrete_time_unit_normed);
    }
    if( is_calculate_cdf ){
      // convert nucsub_cdf from unnormalized distro to cdf:
      for(unsigned i(0);i<NUC3MER_SIZE;++i){
        pdistro_norm(nucsub_cdf_cat+i*(NUC::SIZE-1),nucsub_cdf_cat+(i+1)*(NUC::SIZE-1));
        pdistro_to_cdf_inplace(nucsub_cdf_cat+i*(NUC::SIZE-1),(NUC::SIZE-1));
      }
    }
  }
}




struct ancestral_seq_callinfo {
  const root_gtor& root;
  const unsigned size;
  const unsigned n_cats;
  const std::vector<unsigned>& cat_seq;
};




// group_seq specifies continuous sequence ranges: sequence
// dependencies are reset at the boundary between each range
//
static
void
get_ancestral_seq_nsc4(const ancestral_seq_callinfo& asc,
                       const std::vector<unsigned>& group_seq,
                       const RATE_GTOR_MODEL::index_t c4_type,
                       std::vector<unsigned>& seq){

  struct nsc4_root_cdf {
    nsc4_root_cdf(){
      for(unsigned i(0);i<NUC::SIZE;++i){
        nscodon_cdf_by_nx[i] = nscodon_cdf_by_nx_base+NSC4::get_nscodon_offset(static_cast<NUC::index_t>(i));
      }
    }

    prob_t cdf[NSC4::SIZE];
    prob_t* nscodon_cdf_by_nx[NUC::SIZE];
    prob_t nscodon_cdf_by_nx_base[NSC4::SIZE]; ///< provides storage for nscodon_cdf_by_nx;
  };

  // prepare root_distro pdf's
  nsc4_root_cdf* root_cdf(new nsc4_root_cdf[asc.n_cats]);

  for(unsigned i(0);i<asc.n_cats;++i){
    const prob_t* root_pdistro(asc.root.cat_state_pdistro(i));

    // cdf used for independent 4mer draws:
    pdistro_to_cdf(root_pdistro,root_cdf[i].cdf,NSC4::SIZE);

    // cdf used for draws dependent on the 4mer's extra nuc:
    for(unsigned j(0);j<NSC4::SIZE;++j){ root_cdf[i].nscodon_cdf_by_nx_base[j] = root_pdistro[j]; }

    for(unsigned j(0);j<NUC::SIZE;++j){
      pdistro_norm(root_cdf[i].nscodon_cdf_by_nx[j],root_cdf[i].nscodon_cdf_by_nx[j]+NSCODON::SIZE);
      pdistro_to_cdf_inplace(root_cdf[i].nscodon_cdf_by_nx[j],NSCODON::SIZE);
    }
  }


  seq.resize(asc.size);

  {
    // simulation values for C4PRE:
    int seq_start(0);
    int seq_end(asc.size);
    int seq_increment(1);
    unsigned dependent_nuc(2);

    // simulation values for C4POST:
    if(c4_type == RATE_GTOR_MODEL::C4POST){
      seq_start = asc.size-1;
      seq_end = -1;
      seq_increment = -1;
      dependent_nuc = 0;
    }

    NUC::index_t n[CODON::BASE_SIZE];
    for(int i(seq_start); i != seq_end; i += seq_increment){

      const nsc4_root_cdf& r(root_cdf[asc.cat_seq[i]]);

      if(i==seq_start || group_seq[i-seq_increment] != group_seq[i]){
        // for C4PRE(C4POST): first (last) position in each group is
        // not dependent on nuc from the preceding (following) codon:
        //
        seq[i] = NSC4::decode_nscodon(random_cdf_variate(r.cdf,NSC4::SIZE));
      } else {
        // for C4PRE(C4POST): all other positions are found as a
        // function of the dependent_nuc in the preceding (following)
        // codon:
        //
        NSCODON::decode(n,static_cast<NSCODON::index_t>(seq[i-seq_increment]));
        seq[i] = random_cdf_variate(r.nscodon_cdf_by_nx[n[dependent_nuc]],NSCODON::SIZE);
      }
    }
  }

  delete [] root_cdf;
}




static
void
get_ancestral_seq_nsc5_prepost(const ancestral_seq_callinfo& asc,
                               const std::vector<unsigned>& group_seq,
                               const RATE_GTOR_MODEL::index_t model,
                               std::vector<unsigned>& seq){

  struct nsc5pp_root_cdf {
    nsc5pp_root_cdf(){
      for(unsigned i(0);i<NUC::SIZE;++i){
        for(unsigned j(0);j<NUC::SIZE;++j){
          nscodon_cdf_by_n2x[j+i*NUC::SIZE] =
            nscodon_cdf_by_n2x_base+NSC5::get_nscodon_offset(static_cast<NUC::index_t>(i),static_cast<NUC::index_t>(j));
        }
      }
    }

    prob_t cdf[NSC5::SIZE];
    prob_t* nscodon_cdf_by_n2x[NUC::SIZE*NUC::SIZE];
    prob_t nscodon_cdf_by_n2x_base[NSC5::SIZE]; ///< provides storage for nscodon_cdf_by_nx;
  };

  // prepare root_distro pdf's
  nsc5pp_root_cdf* root_cdf(new nsc5pp_root_cdf[asc.n_cats]);

  for(unsigned i(0);i<asc.n_cats;++i){
    const prob_t* root_pdistro(asc.root.cat_state_pdistro(i));

    // cdf used for independent 5mer draws:
    pdistro_to_cdf(root_pdistro,root_cdf[i].cdf,NSC5::SIZE);

    // cdf used for draws dependent on the 4mer's extra nuc:
    for(unsigned j(0);j<NSC5::SIZE;++j){ root_cdf[i].nscodon_cdf_by_n2x_base[j] = root_pdistro[j]; }

    for(unsigned j(0);j<NUC::SIZE*NUC::SIZE;++j){
      pdistro_norm(root_cdf[i].nscodon_cdf_by_n2x[j],root_cdf[i].nscodon_cdf_by_n2x[j]+NSCODON::SIZE);
      pdistro_to_cdf_inplace(root_cdf[i].nscodon_cdf_by_n2x[j],NSCODON::SIZE);
    }
  }

  seq.resize(asc.size);

  {
    // simulation values for C5PRE:
    int seq_start(0);
    int seq_end(asc.size);
    int seq_increment(1);
    unsigned n1_cpos(1);
    unsigned n2_cpos(2);

    // simulation values for C5POST:
    if(model == RATE_GTOR_MODEL::C5POST){
      seq_start = asc.size-1;
      seq_end = -1;
      seq_increment = -1;
      n1_cpos = 0;
      n2_cpos = 1;
    }

    NUC::index_t n[CODON::BASE_SIZE];
    for(int i(seq_start); i != seq_end; i += seq_increment){

      const nsc5pp_root_cdf& r(root_cdf[asc.cat_seq[i]]);

      if(i==seq_start || group_seq[i-seq_increment] != group_seq[i]){
        // for C4PRE(C4POST): first (last) position in each group is
        // not dependent on nuc from the preceding (following) codon:
        //
        seq[i] = NSC5::decode_nscodon(random_cdf_variate(r.cdf,NSC5::SIZE));
      } else {
        // for C4PRE(C4POST): all other positions are found as a
        // function of the dependent_nuc in the preceding (following)
        // codon:
        //
        NSCODON::decode(n,static_cast<NSCODON::index_t>(seq[i-seq_increment]));
        seq[i] = random_cdf_variate(r.nscodon_cdf_by_n2x[n[n2_cpos]+NUC::SIZE*n[n1_cpos]],NSCODON::SIZE);
      }
    }
  }

  delete [] root_cdf;
}




// group_seq specifies continuous sequence ranges: sequence
// dependencies are reset at the boundary between each range
//
static
void
get_ancestral_seq_nsc5(const ancestral_seq_callinfo& asc,
                       const std::vector<unsigned>& group_seq,
                       std::vector<unsigned>& seq){

  struct nsc5_root_cdf {
    nsc5_root_cdf() {
      for(unsigned i(0);i<DINUC::SIZE;++i){
        for(unsigned j(0);j<TRINUC::SIZE;++j) trinuc_cdf_by_dinuc[i][j] = 0.;
      }
    }

    prob_t cdf[NSC5::SIZE]; ///< cdf used for independent 5mer draws
    prob_t trinuc_cdf_by_dinuc[DINUC::SIZE][TRINUC::SIZE]; ///< cdf used for draws dependent on pre_nuc and n0
  };


  // prepare root_distro pdf's
  nsc5_root_cdf* root_cdf(new nsc5_root_cdf[asc.n_cats]);

  for(unsigned i(0);i<asc.n_cats;++i){
    const prob_t* root_pdistro(asc.root.cat_state_pdistro(i));

    pdistro_to_cdf(root_pdistro,root_cdf[i].cdf,NSC5::SIZE);

    NUC::index_t n[CODON::BASE_SIZE];
    for(unsigned j(0);j<NSC5::SIZE;++j) {
      NUC::index_t n_pre,n_post;
      NSCODON::index_t c;
      NSC5::decode(n_pre,n_post,c,static_cast<NSC5::index_t>(j));
      NSCODON::decode(n,c);
      const DINUC::index_t d = DINUC::encode(n_pre,n[0]);
      const TRINUC::index_t t = TRINUC::encode(n[1],n[2],n_post);

      root_cdf[i].trinuc_cdf_by_dinuc[d][t] = root_pdistro[j];
    }

    for(unsigned j(0);j<DINUC::SIZE;++j){
      pdistro_norm(root_cdf[i].trinuc_cdf_by_dinuc[j],root_cdf[i].trinuc_cdf_by_dinuc[j]+TRINUC::SIZE);
      pdistro_to_cdf_inplace(root_cdf[i].trinuc_cdf_by_dinuc[j],TRINUC::SIZE);
    }
  }


  seq.resize(asc.size);

  NUC::index_t n[CODON::BASE_SIZE+1];
  // initialize to invalid value for safety:
  for(unsigned i(0);i<(CODON::BASE_SIZE+1);++i) n[i] = NUC::N;

  for(int i(0);i<static_cast<int>(asc.size);++i){

    const nsc5_root_cdf& r(root_cdf[asc.cat_seq[i]]);

    if(i==0 || group_seq[i-1] != group_seq[i]){
      // first position in each group is not dependent on nuc from
      // the preceding codon:
      //
      const NSC5::index_t f(random_cdf_variate(r.cdf,NSC5::SIZE));
      const NSCODON::index_t c(NSC5::decode_nscodon(f));
      seq[i] = c;

      NSCODON::decode(n,c);
      n[3] = NSC5::decode_nuc2(f);
    } else {
      const unsigned t = random_cdf_variate(r.trinuc_cdf_by_dinuc[DINUC::encode(n[2],n[3])],TRINUC::SIZE);
      n[0] = n[3];
      TRINUC::decode(n+1,static_cast<TRINUC::index_t>(t));

      seq[i] = NSCODON::encode(n);
    }
  }

  delete [] root_cdf;
}




// group_seq specifies continuous sequence ranges: sequence
// dependencies are reset at the boundary between each range
//
static
void
get_ancestral_seq_dinuc(const ancestral_seq_callinfo& asc,
                        const std::vector<unsigned>& group_seq,
                        std::vector<unsigned>& seq){

  struct dinuc_root_cdf {
    dinuc_root_cdf(){
      for(unsigned i(0);i<NUC::SIZE;++i){
        nuc_cdf_by_nx[i] = cdf+(i*NUC::SIZE);
      }
    }

    prob_t cdf[DINUC::SIZE];
    prob_t* nuc_cdf_by_nx[NUC::SIZE];
    prob_t nuc_cdf_by_nx_base[DINUC::SIZE]; ///< provides storage for nuc_cdf_by_nx;
  };

  // prepare root_distro pdf's
  dinuc_root_cdf* root_cdf(new dinuc_root_cdf[asc.n_cats]);

  for(unsigned i(0);i<asc.n_cats;++i){
    const prob_t* root_pdistro(asc.root.cat_state_pdistro(i));

    pdistro_to_cdf(root_pdistro,root_cdf[i].cdf,DINUC::SIZE);

    for(unsigned j(0);j<DINUC::SIZE;++j){ root_cdf[i].nuc_cdf_by_nx_base[j] = root_pdistro[j]; }

    for(unsigned j(0);j<NUC::SIZE;++j){
      pdistro_norm(root_cdf[i].nuc_cdf_by_nx[j],root_cdf[i].nuc_cdf_by_nx[j]+NUC::SIZE);
      pdistro_to_cdf_inplace(root_cdf[i].nuc_cdf_by_nx[j],NUC::SIZE);
    }
  }


  seq.resize(asc.size);

  {
    const int seq_start(0);
    const int seq_end(asc.size);
    const int seq_increment(1);

    for(int i(seq_start); i != seq_end; i += seq_increment){

      const dinuc_root_cdf& r(root_cdf[asc.cat_seq[i]]);

      if(i==seq_start || group_seq[i-seq_increment] != group_seq[i]){
        // first position in each group is not dependent on nuc from
        // the preceding dinuc:
        //
        seq[i] = DINUC::decode_n0(random_cdf_variate(r.cdf,DINUC::SIZE));
      } else {
        // all other positions are found as a function of the
        // dependent_nuc in the preceding dinuc:
        //
        int nx = DINUC::decode_nx(seq[i-seq_increment]);
        seq[i] = random_cdf_variate(r.nuc_cdf_by_nx[nx],NUC::SIZE);
      }
    }
  }

  delete [] root_cdf;
}




// group_seq specifies continuous sequence ranges: sequence
// dependencies are reset at the boundary between each range
//
static
void
get_ancestral_seq_trinuc(const ancestral_seq_callinfo& asc,
                         const std::vector<unsigned>& group_seq,
                         std::vector<unsigned>& seq){

  struct trinuc_root_cdf {
    prob_t cdf[TRINUC::SIZE];
    prob_t dinuc_cdf_by_nuc0[NUC::SIZE][DINUC::SIZE];
    prob_t nuc_cdf_by_dinuc0[DINUC::SIZE][NUC::SIZE];
  };


  // prepare root_distro pdf's
  trinuc_root_cdf* root_cdf(new trinuc_root_cdf[asc.n_cats]);

  for(unsigned i(0);i<asc.n_cats;++i){
    const prob_t* root_pdistro(asc.root.cat_state_pdistro(i));

    pdistro_to_cdf(root_pdistro,root_cdf[i].cdf,TRINUC::SIZE);

    NUC::index_t nuc[TRINUC::BASE_SIZE];
    for(unsigned j(0);j<TRINUC::SIZE;++j){
      TRINUC::decode(nuc,static_cast<TRINUC::index_t>(j));
      root_cdf[i].dinuc_cdf_by_nuc0[nuc[0]][DINUC::encode(nuc[1],nuc[2])] = root_pdistro[j];
    }
    for(unsigned j(0);j<NUC::SIZE;++j){
      pdistro_norm(root_cdf[i].dinuc_cdf_by_nuc0[j],root_cdf[i].dinuc_cdf_by_nuc0[j]+DINUC::SIZE);
      pdistro_to_cdf_inplace(root_cdf[i].dinuc_cdf_by_nuc0[j],DINUC::SIZE);
    }

    for(unsigned j(0);j<TRINUC::SIZE;++j){
      TRINUC::decode(nuc,static_cast<TRINUC::index_t>(j));
      root_cdf[i].nuc_cdf_by_dinuc0[DINUC::encode(nuc[0],nuc[1])][nuc[2]] = root_pdistro[j];
    }
    for(unsigned j(0);j<DINUC::SIZE;++j){
      pdistro_norm(root_cdf[i].nuc_cdf_by_dinuc0[j],root_cdf[i].nuc_cdf_by_dinuc0[j]+NUC::SIZE);
      pdistro_to_cdf_inplace(root_cdf[i].nuc_cdf_by_dinuc0[j],NUC::SIZE);
    }
  }


  seq.resize(asc.size);

  {
    const int seq_start(0);
    const int seq_end(asc.size);
    const int seq_increment(1);

    NUC::index_t nuc[TRINUC::SIZE];
    for(int i(seq_start); i != seq_end; i += seq_increment){

      const trinuc_root_cdf& r(root_cdf[asc.cat_seq[i]]);

      if(i==seq_start || group_seq[i-seq_increment] != group_seq[i]){
        // first position in each group is not dependent on nuc from
        // the preceding dinuc:
        //
        TRINUC::decode(nuc,static_cast<TRINUC::index_t>(random_cdf_variate(r.cdf,TRINUC::SIZE)));
        seq[i] = nuc[2];
      } else if((i-1)==seq_start || group_seq[(i-1)-seq_increment] != group_seq[i-1]) {
        //
        //
        seq[i] = DINUC::decode_n0(random_cdf_variate(r.dinuc_cdf_by_nuc0[seq[i-seq_increment]],DINUC::SIZE));
      } else {
        // all other positions are found as a function of the
        // dependent_dinuc in the preceding trinuc:
        //
        const NUC::index_t nx2 = static_cast<NUC::index_t>(seq[(i-1)-seq_increment]);
        const NUC::index_t nx1 = static_cast<NUC::index_t>(seq[i-seq_increment]);
        seq[i] = random_cdf_variate(r.nuc_cdf_by_dinuc0[DINUC::encode(nx2,nx1)],NUC::SIZE);
      }
    }
  }

  delete [] root_cdf;
}



const unsigned left_nuc_pad_size(2);
const unsigned right_nuc_pad_size(2);
const unsigned nuc_pad_size(left_nuc_pad_size+right_nuc_pad_size);



static
void
sim_data_nuc_to_nuc_seq_data(const bi_tree& tree,
                             const std::vector<std::vector<unsigned > >& treenode_seq,
                             const unsigned n_groups,
                             const prob_t assigned_cat_prob,
                             const std::vector<unsigned>& cat_seq,
                             const std::vector<unsigned>& group_seq,
                             nuc_seq_data& nsd){

  const unsigned n_leaves(tree.leaf_size());
  group_range gr(group_seq,left_nuc_pad_size,right_nuc_pad_size);

  // chop up present day sequences into nuc format:
  for(unsigned g(0);g<n_groups;++g){
    unsigned group_start,group_stop;
    if(! gr.get_gbounds(g,group_start,group_stop)) continue;

    const unsigned group_size(group_stop-group_start);
    align_dat_t& na(nsd.dat[g]);
    na.seq.init(n_leaves,group_size);
    na.info.resize(group_size);

    for(unsigned i(0);i<group_size;++i){
      for(unsigned t(0);t<n_leaves;++t){
        const std::vector<unsigned>& leaf_seq(treenode_seq[tree.leaf_node(t)->node_id()]);
        na.seq[t][i] = static_cast<NUC::index_t>(leaf_seq[group_start+i]);
      }

      align_col_t& nai(na.info[i]);

      nai.class_no=UNASSIGNED_CAT_ID;
      if(assigned_cat_prob>0.){
        const prob_t r(random_uniform());
        if(r<assigned_cat_prob){ nai.class_no=cat_seq[group_start+i]+1; }
      }

      nai.is_continuous = (i!=0);
    }
  }
}


#ifdef DIRECT_SITE_CONV
static
void
accumulate_simulated_seq_group_nuc(const bi_tree& tree,
                                   const std::vector<std::vector<unsigned > >& treenode_seq,
                                   const RATE_GTOR_MODEL::index_t rgm,
                                   const prob_t assigned_cat_prob,
                                   const std::vector<unsigned>& cat_seq,
                                   const std::vector<unsigned>& group_seq,
                                   site_data& sd){

  assert(RATE_GTOR_MODEL::base_size_conditioned(rgm)==1);

  const unsigned n_leaves(tree.leaf_size());
  const unsigned full_size(treenode_seq[0].size());
  const SITE_MODEL::index_t smodel(RATE_GTOR_MODEL::convert_to_site_model(rgm));
  const unsigned ambig_state(SITE_MODEL::ambig_state(smodel));

  site_code sc(n_leaves);

  // chop up present day sequences into cq data format:
  for(unsigned i(0);i<(full_size);++i){

    // skip nuc left and right pads in every group
    bool is_skip(false);

    if(i<left_nuc_pad_size || (i+right_nuc_pad_size)>=full_size) {
      is_skip = true;
    } else {
      const unsigned this_group(group_seq[i]);
      for(unsigned j(0);j<left_nuc_pad_size;++j){
        if(group_seq[i-(j+1)] != this_group) is_skip=true;
      }
      for(unsigned j(0);j<right_nuc_pad_size;++j){
        if(group_seq[i+(j+1)] != this_group) is_skip=true;
      }
    }

    if(is_skip) continue;


    const int offset(RATE_GTOR_MODEL::base_size_repeat_offset(rgm));
    NUC::index_t nuc[SITE_MODEL::MAX_BASE_SIZE];

    for(unsigned t(0);t<n_leaves;++t){
      const std::vector<unsigned>& leaf_seq(treenode_seq[tree.leaf_node(t)->node_id()]);
      for(int j(0);j<(1-offset);++j){
        nuc[j] = static_cast<NUC::index_t>(leaf_seq[i+offset+j]);
      }
      sc.set_taxid(t,SITE_MODEL::encode_nuc(smodel,nuc));
    }

    unsigned cat_no(UNASSIGNED_CAT_ID);
    if(assigned_cat_prob>0.){
      const prob_t r(random_uniform());
      if(r<assigned_cat_prob){ cat_no=cat_seq[i]+1; }
    }

    const unsigned group_no(group_seq[i]);

    sd.site_count[sc][count_data_code(group_no,cat_no)]++;
  }
}
#endif


// pos = 1-indexed position of change in the codon
//
// this surveys the possible nucleotide changes for one background
// [nuc-threemer+codon=3or4 nucs] at one position. n[1,2,3] are
// interpreted as the codon n[0] is interpreted as the preceding
// nucleotide for pos == 1 or the following nucleotide if pos == 3
//
static
void
codon_mut_prob_context_position_type(const rate_gtor_nscodon_base& r,
                                     const NUC::index_t n[],
                                     prob_t& subprob,
                                     prob_t* subdist,
                                     const unsigned pos,
                                     const unsigned cat,
                                     const unsigned branch_cat_set,
                                     const bool is_include_codon_selection){

  static const unsigned NMER_SIZE(4);
  if(pos == 0) die("invalid pos arg to get_discrete_time_prob_codon_context_position_type");

  const NSCODON::index_t codon = NSCODON::encode(n+1);
  const NSAA::index_t aa = codon_trans_known(codon);
  const NUC::index_t np = n[pos];
  NUC::index_t ncopy[NMER_SIZE];
  for(unsigned i(0);i<NMER_SIZE;++i) ncopy[i] = n[i];

  // loop through possible nuc changes at pos:
  //
  for(unsigned j(0);j<NUC::SIZE;++j){
    const NUC::index_t nj(static_cast<NUC::index_t>(j));
    if(nj == np) continue;
    ncopy[pos] = nj;

    const NSCODON::index_t codon_j(NSCODON::encode(ncopy+1));

    smlfloat this_rate(r.cat_nuc_context_mut_rate(ncopy[pos-1],np,ncopy[(pos+1)%NMER_SIZE],
                                                  nj,cat,branch_cat_set,0,0,codon_j));

    if( is_include_codon_selection ){
      if( codon_j == NSCODON::NNN ){ // we've hit a stop codon
        this_rate = 0.;
      } else {
        const NSAA::index_t aa_j = codon_trans_known(codon_j);
        this_rate *= r.codon_bias(codon,codon_j);
        if(aa != aa_j) this_rate *= r.cat_nsaa_param(aa,aa_j,cat,branch_cat_set);
      }
    }

    subdist[nj+(nj<np?0:-1)] = this_rate;
  }
  subprob = array_sum(subdist,(NUC::SIZE-1));
}



// pos = 1-indexed position of change in the codon
//
static
void
codon_mut_prob_context_position(const rate_gtor_nscodon_base& r,
                                prob_t* codonsub,
                                prob_t* codonsub_cdf,
                                const unsigned pos,
                                const unsigned cat,
                                const unsigned branch_cat_set,
                                const bool is_include_codon_selection,
                                const bool is_calculate_cdf){

  static const unsigned NMER_SIZE(4);
  if(pos<1 || pos>3) die("invalid pos arg to get_discrete_time_prob_codon_context_position");

  for(unsigned i(0);i<NSC4::SIZE;++i){
    NUC::index_t n[NMER_SIZE];
    {
      NSCODON::index_t c;
      NSC4::decode(n[0],c,i);
      NSCODON::decode(n[1],n[2],n[3],c);
    }

    prob_t subprob;
    prob_t subdist[NUC::SIZE-1];
    codon_mut_prob_context_position_type(r,n,subprob,subdist,pos,
                                         cat,branch_cat_set,is_include_codon_selection);

    // load local data structures (subprob and subdist) into exported
    // data structures:
    for(unsigned j(0);j<NUC::SIZE;++j){
      unsigned ns5mer;
      if      ( pos == 1 || pos == 2 ) {
        const NSC4::index_t pre4mer(i);
        ns5mer = pre4mer+j*NSC4::SIZE;
      }else{ // pos == 3
        const NSCODON::index_t c = NSCODON::encode(n[1],n[2],n[3]);
        const NSC4::index_t pre4mer = NSC4::encode(static_cast<NUC::index_t>(j),c);
        ns5mer = pre4mer+n[0]*NSC4::SIZE;
      }

      // add up codonsub over all positions
      codonsub[ns5mer] += subprob;

      if( is_calculate_cdf ){
        prob_t* codonsub_cdf_5mer = codonsub_cdf+ns5mer*CODON_CDF_SIZE;
        prob_t* codonsub_cdf_5mer_pos = codonsub_cdf_5mer+(pos-1)*(NUC::SIZE-1);
        for(unsigned k(0);k<(NUC::SIZE-1);++k){
          codonsub_cdf_5mer_pos[k] = subdist[k];
        }
      }
    }
  }
}




// if is_use_discrete_time, time is supplied in subs per site and
// normed according to the rate matrix
//
// a second run with selection turned off is used as part of
// calculating simulation time
//
static
void
codon_mut_prob_context(const rate_gtor_nscodon_base& r,
                       prob_t** codonsub,       // [n_cats][NS5MER_SIZE]
                       prob_t** codonsub_cdf,   // [n_cats][NS5MER_SIZE*CODON_CDF_SIZE]
                       const bool is_use_discrete_time,
                       const smlfloat discrete_time_unit_normed,
                       const unsigned n_cats,
                       const time_gtor& tg,
                       const unsigned branch_no,
                       const cat_manager& cm,
                       bool is_include_codon_selection = true,
                       const bool is_calculate_cdf = true){

  const smlfloat etime(tg.cat_expected_branch_time(branch_no));

  for(unsigned c(0);c<n_cats;++c) std::fill(codonsub[c],codonsub[c]+NS5MER_SIZE,0.);

  for(unsigned c(0);c<n_cats;++c){
    const unsigned bcs(cm.get_branch_cat_set(c,branch_no));

    prob_t* codonsub_cat(codonsub[c]);
    prob_t* codonsub_cdf_cat(0);
    if( is_calculate_cdf ) codonsub_cdf_cat = codonsub_cdf[c];

    for(unsigned pos(1);pos<=3;++pos){
      codon_mut_prob_context_position(r,codonsub_cat,codonsub_cdf_cat,
                                      pos,c,bcs,is_include_codon_selection,is_calculate_cdf);
    }

    const smlfloat branch_time(tg.branch_time(branch_no,c));
    const smlfloat branch_cat_scale(branch_time <= 0. ? 0. : branch_time/etime);

    array_scale(codonsub_cat,NS5MER_SIZE,branch_cat_scale);

    if(is_use_discrete_time){
      array_scale(codonsub_cat,NS5MER_SIZE,discrete_time_unit_normed);
    }
    if( is_calculate_cdf ){
      // convert codonsub_cdf from unnormalized distro to cdf:
      for(unsigned i(0);i<NS5MER_SIZE;++i){
        pdistro_norm(codonsub_cdf_cat+i*CODON_CDF_SIZE,codonsub_cdf_cat+(i+1)*CODON_CDF_SIZE);
        pdistro_to_cdf_inplace(codonsub_cdf_cat+i*CODON_CDF_SIZE,CODON_CDF_SIZE);
      }
    }
  }
}




const unsigned left_codon_pad_size(1);
const unsigned right_codon_pad_size(1);
const unsigned codon_pad_size(left_codon_pad_size+right_codon_pad_size);


static
void
sim_data_codon_to_nuc_seq_data(const bi_tree& tree,
                               const std::vector<std::vector<unsigned > >& treenode_seq,
                               const unsigned n_groups,
                               const prob_t assigned_cat_prob,
                               const std::vector<unsigned>& cat_seq,
                               const std::vector<unsigned>& group_seq,
                               nuc_seq_data& nsd){

  const unsigned n_leaves(tree.leaf_size());
  group_range gr(group_seq,left_codon_pad_size,right_codon_pad_size);

  // chop up present day sequences into nuc format:
  for(unsigned g(0);g<n_groups;++g){
    unsigned group_start,group_stop;
    if(! gr.get_gbounds(g,group_start,group_stop)) continue;

    const unsigned group_size(group_stop-group_start);
    align_dat_t& na(nsd.dat[g]);
    na.seq.init(n_leaves,group_size*CODON::BASE_SIZE);
    na.info.resize(group_size*CODON::BASE_SIZE);

    for(unsigned i(0);i<group_size;++i){
      for(unsigned t(0);t<n_leaves;++t){
        const std::vector<unsigned>& leaf_seq(treenode_seq[tree.leaf_node(t)->node_id()]);
        NSCODON::decode(na.seq[t]+i*CODON::BASE_SIZE,
                        static_cast<NSCODON::index_t>(leaf_seq[group_start+i]));
      }

      unsigned codon_cat_no(UNASSIGNED_CAT_ID);
      if(assigned_cat_prob>0.){
        const prob_t r(random_uniform());
        if(r<assigned_cat_prob) codon_cat_no=cat_seq[group_start+i]+1;
      }

      for(unsigned j(0);j<CODON::BASE_SIZE;++j){
        const unsigned nuci((i*CODON::BASE_SIZE)+j);
        align_col_t& nai(na.info[nuci]);

        nai.class_no=codon_cat_no;
        nai.codon_pos=j;
        nai.is_continuous = (nuci!=0);
      }
    }
  }
}


#ifdef DIRECT_SITE_CONV
static
void
accumulate_simulated_seq_group(const bi_tree& tree,
                               const std::vector<std::vector<unsigned > >& treenode_seq,
                               const RATE_GTOR_MODEL::index_t rgm,
                               const prob_t assigned_cat_prob,
                               const std::vector<unsigned>& cat_seq,
                               const std::vector<unsigned>& group_seq,
                               site_data& sd){

  assert(RATE_GTOR_MODEL::base_size_conditioned(rgm)==CODON::BASE_SIZE);

  const unsigned n_leaves(tree.leaf_size());
  const unsigned full_size(treenode_seq[0].size());
  const SITE_MODEL::index_t smodel(RATE_GTOR_MODEL::convert_to_site_model(rgm));
  const unsigned ambig_state(SITE_MODEL::ambig_state(smodel));

  site_code sc(n_leaves);

  // chop up present day sequences into cq data format:
  for(unsigned i(0);i<full_size;++i){

    // skip nuc left and right pads in every group
    bool is_skip(false);

    if(i<left_codon_pad_size || (i+right_codon_pad_size)>=full_size) {
      is_skip = true;
    } else {
      const unsigned this_group(group_seq[i]);
      for(unsigned j(0);j<left_codon_pad_size;++j){
        if(group_seq[i-(j+1)] != this_group) is_skip=true;
      }
      for(unsigned j(0);j<right_codon_pad_size;++j){
        if(group_seq[i+(j+1)] != this_group) is_skip=true;
      }
    }

    if(is_skip) continue;

    const int offset(RATE_GTOR_MODEL::base_size_repeat_offset(rgm));

    NUC::index_t nuc[CODON::BASE_SIZE*(1+codon_pad_size)];
    const NUC::index_t* nuc_start(nuc+CODON::BASE_SIZE*left_codon_pad_size);

    for(unsigned t(0);t<n_leaves;++t){
      const std::vector<unsigned>& leaf_seq(treenode_seq[tree.leaf_node(t)->node_id()]);
      for(unsigned j(0);j<(1+codon_pad_size);++j){
        const int codon_offset(j-left_codon_pad_size);
        NSCODON::decode(nuc+CODON::BASE_SIZE*(j),static_cast<NSCODON::index_t>(leaf_seq[i+codon_offset]));
      }
      sc.set_taxid(t,SITE_MODEL::encode_nuc(smodel,nuc_start+offset));
    }

    unsigned cat_no(UNASSIGNED_CAT_ID);

    if(assigned_cat_prob>0.){
      const prob_t r(random_uniform());
      if(r<assigned_cat_prob){ cat_no=cat_seq[i]+1; }
    }
    const unsigned group_no(group_seq[i]);

    sd.site_count[sc][count_data_code(group_no,cat_no)]++;
  }
}
#endif



void
simulate_data_context(const sim_options& sim_opt,
                      const subs_ml_model& mdl,
                      nuc_seq_data& nsd){

  const cat_manager& cm(mdl.get_cat_manager());
  const unsigned n_cats(cm.cat_size());

  const rate_gtor_nuc_base* rgn_ptr(dynamic_cast<const rate_gtor_nuc_base*>(&mdl.get_rate_gtor()));
  if(rgn_ptr==0){
    pass_away("simulate_data_context requires a nucleotide-based model");
  }
  const rate_gtor_nuc_base& rgn(*rgn_ptr);

  const rate_gtor_nscodon_base* rgc_ptr(dynamic_cast<const rate_gtor_nscodon_base*>(&rgn));

  const bool is_codon_site_model(rgc_ptr != 0);

  // site type of simulated model
  const RATE_GTOR_MODEL::index_t rgm(rgn.model_type());

  //nuc site default:
  unsigned context_state_size(NUC3MER_SIZE);
  unsigned mutation_state_size(NUC::SIZE-1);
  unsigned time_nuc_factor(1);
  if(is_codon_site_model){
    context_state_size = NS5MER_SIZE;
    mutation_state_size = CODON_CDF_SIZE;
    time_nuc_factor=CODON::BASE_SIZE;
  }

  // use background nuc distro to handle the edges of each group:
  simple_matrix<prob_t> bg_nuc_cdf_cat(n_cats,NUC::SIZE);
  for(unsigned c(0);c<n_cats;++c){
    pdistro_to_cdf(rgn.bg_pdistro_nuc_cat(c),bg_nuc_cdf_cat[c],NUC::SIZE);
  }

  // assign sites to groups:
  std::vector<unsigned> group_seq;
  get_sim_group_seq(group_seq,sim_opt);

  const unsigned n_groups(group_seq.back()+1);

  // 3. if simulating 4mers/dinuc(5mers/trinuc), add one (two) extra codon(s) for each group
  unsigned full_size(0);
  std::vector<unsigned> full_group_seq;
  {
    // set a uniform pad for all codon site types or all nuc site types, so that sim can be repeated
    // with the same seed:
    //
    const unsigned pad_size(is_codon_site_model ? codon_pad_size : nuc_pad_size );

    full_size = sim_opt.size+(pad_size*n_groups);
    for(unsigned i(0);i<sim_opt.size;++i) {
      unsigned copy_num(1);
      if(i==0 || group_seq[i-1] != group_seq[i]) copy_num=1+pad_size;
      for(unsigned k(0);k<copy_num;++k) {
        full_group_seq.push_back(group_seq[i]);
      }
    }
  }

  // assigned cats to sites
  std::vector<unsigned> cat_seq;
  get_sim_cat_seq(cat_seq,full_group_seq,full_size,mdl);


  // get_ancestral seq:
  const bi_tree& tree(mdl.tree());
  const unsigned n_nodes(tree.node_size());
  const unsigned n_branches(tree.branch_size());
  std::vector<std::vector<unsigned> > treenode_seq(n_nodes);

  {
    const unsigned root_id(tree.root()->node_id());
    std::vector<unsigned>& root_seq(treenode_seq[root_id]);
    const ancestral_seq_callinfo asc = { mdl.get_root_gtor(), full_size, n_cats, cat_seq };

    using namespace RATE_GTOR_MODEL;
    switch(rgm){
    case C4PRE:
    case C4POST:
      get_ancestral_seq_nsc4(asc,full_group_seq,rgm,root_seq);
      break;
    case RATE_GTOR_MODEL::DINUC:
      get_ancestral_seq_dinuc(asc,full_group_seq,root_seq);
      break;
    case RATE_GTOR_MODEL::TRINUC:
      get_ancestral_seq_trinuc(asc,full_group_seq,root_seq);
      break;
    case C5PRE:
    case C5POST:
      get_ancestral_seq_nsc5_prepost(asc,full_group_seq,rgm,root_seq);
      break;
    case C5:
      get_ancestral_seq_nsc5(asc,full_group_seq,root_seq);
      break;
    default:
      get_ancestral_seq(mdl,full_size,cat_seq,root_seq);
      break;
    }
  }

  // sim time data structures:
  simple_init_matrix3d<prob_t> sitesub_prob_neutral;
  std::vector<std::vector<smlfloat> > cat_sim_time(n_branches);
  for(unsigned i(0);i<n_branches;++i){ cat_sim_time[i].resize(n_cats); }

  if(sim_opt.is_report_time){
    if(!is_codon_site_model) die("sim time routine not implemented for noncoding models");

    sitesub_prob_neutral.init(n_branches,n_cats,NS5MER_SIZE);
  }

  // primary simulation probabilities
  simple_matrix3d<prob_t> sitesub_prob(n_branches,n_cats,context_state_size);
  simple_matrix3d<prob_t> sitesub_cdf(n_branches,n_cats,context_state_size*mutation_state_size);

  if       (sim_opt.method == SIM_MODEL::CONTINUOUS){
    if(!is_codon_site_model) die("continuous simulator written for coding sequence models only");

    if(sim_opt.is_report_time){
      for(unsigned b(0);b<n_branches;++b){
        codon_mut_prob_context(*rgc_ptr,sitesub_prob_neutral[b],0,false,0,
                               n_cats,mdl.get_time_gtor(),b,cm,false,false);
      }
    }

    for(unsigned b(0);b<n_branches;++b){
      codon_mut_prob_context(*rgc_ptr,sitesub_prob[b],sitesub_cdf[b],false,0,
                             n_cats,mdl.get_time_gtor(),b,cm);
    }

    // simulate
    for(unsigned b(0);b<n_branches;++b){
      const bi_tree_node* branch_node(tree.branch_node(b));
      const smlfloat etime(mdl.get_time_gtor().cat_expected_branch_time(b));
      const prob_t* const * spn_neutral(sitesub_prob_neutral.dim1() == 0 ? 0 : sitesub_prob_neutral[b]);
      simulate_continuous_time_branch_context(treenode_seq[branch_node->parent()->node_id()],cat_seq,
                                              full_group_seq,treenode_seq[branch_node->node_id()],
                                              etime,n_cats,
                                              sitesub_prob[b],
                                              sitesub_cdf[b],bg_nuc_cdf_cat.ptr(),
                                              cat_sim_time[b],
                                              is_codon_site_model,
                                              sim_opt.is_report_time,
                                              spn_neutral);
    }

  } else if(sim_opt.method == SIM_MODEL::DISCRETE) {
    // get subs prob structures for each discrete time unit
    //
    // assuming only one substitution can occur per codon, get the
    // probability of any substitution happening conditioned on the
    // nucleotide 5-mer (the codon +/- 1 nuc in each direction) --
    // then make two separate prob structures, one is the cdf
    // (sitesub_cdf) for substitution type, given that a substitution
    // will happen within the inner codon of the 5-mer. the second is
    // the probability of any change occurring at all within this
    // codon(sitesub_prob).
    //
    const smlfloat time_unit_normed(SIMULATE_DISCRETE_TIME_UNIT/
                                    static_cast<smlfloat>(time_nuc_factor));

    if(sim_opt.is_report_time){
      for(unsigned b(0);b<n_branches;++b){
        codon_mut_prob_context(*rgc_ptr,sitesub_prob_neutral[b],0,true,time_unit_normed,
                               n_cats,mdl.get_time_gtor(),b,cm,false,false);
      }
    }

    for(unsigned b(0);b<n_branches;++b){
      if(is_codon_site_model){
        codon_mut_prob_context(*rgc_ptr,sitesub_prob[b],sitesub_cdf[b],true,time_unit_normed,
                               n_cats,mdl.get_time_gtor(),b,cm);
      } else {
        nuc_mut_prob_context(rgn,sitesub_prob[b],sitesub_cdf[b],true,time_unit_normed,
                             n_cats,mdl.get_time_gtor(),b,cm);
      }
    }

    for(unsigned b(0);b<n_branches;++b){
      const bi_tree_node* branch_node(tree.branch_node(b));
      const smlfloat etime(mdl.get_time_gtor().cat_expected_branch_time(b));
      const prob_t* const * spn_neutral(sitesub_prob_neutral.dim1() == 0 ? 0 : sitesub_prob_neutral[b]);
      simulate_discrete_time_branch_context(treenode_seq[branch_node->parent()->node_id()],
                                            cat_seq,
                                            full_group_seq,treenode_seq[branch_node->node_id()],
                                            etime,
                                            time_unit_normed,
                                            sitesub_prob[b],
                                            sitesub_cdf[b],
                                            bg_nuc_cdf_cat.ptr(),
                                            cat_sim_time[b],
                                            is_codon_site_model,
                                            sim_opt.is_report_time,
                                            spn_neutral);
    }
  }

  // clear memory for internal node sequences, which we won't be using:
  //
  /// \todo simulation storage can be optimized down to
  /// seq_len*n_leaves if storage becomes a bigger issue
  //
  for(unsigned n(0);n<n_nodes;++n){
    if(! tree.node(n)->is_leaf()) treenode_seq[n].clear();
  }

  if(sim_opt.is_report_time) report_sim_tree_times(cm,tree,cat_sim_time,log_os);

#ifdef DIRECT_SITE_CONV
  // site type of simulator output
  const RATE_GTOR_MODEL::index_t output_rgm(sim_opt.output_site_model);

  // write sim data out to sd:
  //
  sim_site_data_init(tree,output_rgm,sd);

  static const char group_label[] = "sim_context_group_";
  for(unsigned i(0);i<n_groups;++i){
    std::ostringstream oss;
    oss << group_label << i;
    sd.group_label.push_back(oss.str());
  }

  // 4. chop simulation data down to (site) size::
  if(is_codon_site_model){
    accumulate_simulated_seq_group(tree,treenode_seq,output_rgm,sim_opt.assigned_cat_prob,
                                   cat_seq,full_group_seq,sd);
  } else {
    accumulate_simulated_seq_group_nuc(tree,treenode_seq,output_rgm,sim_opt.assigned_cat_prob,
                                       cat_seq,full_group_seq,sd);
  }
#else
  // write sim data out to nuc_seq_data:
  //
  sim_seq_data_init(tree,n_groups,n_cats,"sim_context",cm,nsd);

  if(is_codon_site_model){
    sim_data_codon_to_nuc_seq_data(tree,treenode_seq,n_groups,
                                   sim_opt.assigned_cat_prob,cat_seq,full_group_seq,nsd);
  } else {
    sim_data_nuc_to_nuc_seq_data(tree,treenode_seq,n_groups,sim_opt.
                                 assigned_cat_prob,cat_seq,full_group_seq,nsd);
  }
#endif
}

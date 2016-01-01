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

// $Id: lhood_site_cached_estep_baz.cc 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "lhood_site_cached_estep_baz.h"
#include "util/general/die.h"

#ifdef INCLUDE_BAZ_HACK

#include "bioseq_report.h"
#include "lhood_site_cached_estep_shared.h"
#include "nsaa_maps.h"
#include "util/bio/bioseq_util.h"
#include "util/bio/bioseq_util_pdf_conversions.h"
#include "util/math/matrix_util_io.h"

#include <ostream>


struct baz_prob {
  baz_prob() : two_on_leaf1(0.), two_on_leaf2(0.), one_on_each(0.) {}

  void
  increment(prob_t* int_codon_pdistro,
            const NSCODON::index_t leaf1_codon,
            const NSCODON::index_t leaf2_codon,
            const NSCODON::index_t int_codon1,
            const NSCODON::index_t int_codon2){

    two_on_leaf2 += int_codon_pdistro[leaf1_codon];
    two_on_leaf1 += int_codon_pdistro[leaf2_codon];

    one_on_each += int_codon_pdistro[int_codon1]+int_codon_pdistro[int_codon2];
  }

  void
  report(const char* global_label,
         const char* leaf1_label,
         const char* leaf2_label,
         const char* each_label,
         std::ostream& os) const {

    os << "leaf1 2ns_" << global_label << " " << leaf1_label << " " << two_on_leaf1 << "\n";
    os << "leaf2 2ns_" << global_label << " " << leaf2_label << " " << two_on_leaf2 << "\n";
    os << "each_ 2ns_" << global_label << " " << each_label << " " << one_on_each << "\n";
    os << "\n";
  }

  prob_t two_on_leaf1;
  prob_t two_on_leaf2;
  prob_t one_on_each;
};

enum bp_types {
  BP_BASIC, BP_13, BP_AR, BP_AI,
  BP_AN,    BP_HP, BP_HH, BP_PP,
  BP_SIZE
};


struct baz_counts {
  baz_counts()
    : bp(new baz_prob[BP_SIZE]),

      Kn_leaf1_highpath_1(0.),
      Kn_leaf1_highpath_2(0.),
      Kn_leaf2_highpath_1(0.),
      Kn_leaf2_highpath_2(0.),

      Kn_leaf1_highpath_top(0.),
      Kn_leaf1_highpath_bot(0.),
      Kn_leaf2_highpath_top(0.),
      Kn_leaf2_highpath_bot(0.),

      Kn_leaf1_highpath_diff(0.),
      Kn_leaf2_highpath_diff(0.),

      Kn_leaf1_highpath_top_frac(0.),
      Kn_leaf2_highpath_top_frac(0.),

      Kn_leaf1_lowpath_1(0.),
      Kn_leaf1_lowpath_2(0.),
      Kn_leaf2_lowpath_1(0.),
      Kn_leaf2_lowpath_2(0.),

      Kn_leaf1_1(0.),
      Kn_leaf1_2(0.),
      Kn_leaf2_1(0.),
      Kn_leaf2_2(0.),

      Kn_one_on_1(0.),
      Kn_one_on_2(0.),

      two_on_leaf1_parsi(0),
      two_on_leaf2_parsi(0),
      one_on_each_parsi(0),
      one_on_leaf1(0.),
      one_on_leaf2(0.),
      one_syn_on_leaf1(0.),
      one_syn_on_leaf2(0.) {

    for(unsigned i(0);i<(NSAA::SIZE*NSAA::SIZE);++i){
      Kn_highpath_top_sub[i] = 0.;
      Kn_highpath_bot_sub[i] = 0.;
      Kn_highpath_both_sub[i] = 0.;
      Kn_single_sub[i] = 0.;

      two_on_leaf_XX[i] = 0.;
      one_on_each_XX[i] = 0.;
    }

    Kn_highpath_both_sub_ARhigh = 0.;
    Kn_highpath_both_sub_ARlow = 0.;
  }


  ~baz_counts() { delete [] bp; }

  baz_prob* bp;

  smlfloat Kn_leaf1_highpath_1;
  smlfloat Kn_leaf1_highpath_2;
  smlfloat Kn_leaf2_highpath_1;
  smlfloat Kn_leaf2_highpath_2;

  smlfloat Kn_leaf1_highpath_top;
  smlfloat Kn_leaf1_highpath_bot;
  smlfloat Kn_leaf2_highpath_top;
  smlfloat Kn_leaf2_highpath_bot;

  smlfloat Kn_highpath_top_sub[NSAA::SIZE*NSAA::SIZE];
  smlfloat Kn_highpath_bot_sub[NSAA::SIZE*NSAA::SIZE];
  smlfloat Kn_highpath_both_sub[NSAA::SIZE*NSAA::SIZE];
  smlfloat Kn_single_sub[NSAA::SIZE*NSAA::SIZE];

  smlfloat Kn_highpath_both_sub_ARhigh;
  smlfloat Kn_highpath_both_sub_ARlow;

  smlfloat Kn_leaf1_highpath_diff;
  smlfloat Kn_leaf2_highpath_diff;

  smlfloat Kn_leaf1_highpath_top_frac;
  smlfloat Kn_leaf2_highpath_top_frac;

  smlfloat Kn_leaf1_lowpath_1;
  smlfloat Kn_leaf1_lowpath_2;
  smlfloat Kn_leaf2_lowpath_1;
  smlfloat Kn_leaf2_lowpath_2;

  smlfloat Kn_leaf1_1;
  smlfloat Kn_leaf1_2;
  smlfloat Kn_leaf2_1;
  smlfloat Kn_leaf2_2;

  smlfloat Kn_one_on_1;
  smlfloat Kn_one_on_2;

  smlfloat two_on_leaf_XX[NSAA::SIZE*NSAA::SIZE];
  smlfloat one_on_each_XX[NSAA::SIZE*NSAA::SIZE];

  unsigned two_on_leaf1_parsi;
  unsigned two_on_leaf2_parsi;
  unsigned one_on_each_parsi;

  prob_t one_on_leaf1;
  prob_t one_on_leaf2;

  prob_t one_syn_on_leaf1;
  prob_t one_syn_on_leaf2;
};


static
bool
is_hphobic(const NSAA::index_t aa){

  using namespace NSAA;

  switch(aa){
  case G :
  case A :
  case V :
  case L :
  case I :
  case M :
  case Y :
  case F :
  case W :
    return true;
  default :
    return false;
  }
}


static
bool
is_polar(const NSAA::index_t aa){
  return ! is_hphobic(aa);
}


struct tree_walk_cache_estep_baz;


static
void
local_increment_site_stats(const tree_lhood_info& ti,
                           const site_data_fastlup& sdf,
                           const tree_walk_cache_estep_baz& twe,
                           const int site_id,
                           const site_data_fastlup_core::index_type* seq_state,
                           baz_counts& bc);


// info gathered from all org levels. for each site record the cache for the branch above each node
// in both the up and down direction, for all nodes where this info is required for the E-step:
struct tree_walk_cache_estep_baz : public tree_walk_cache_estep_base {

  typedef tree_walk_cache_estep_base base_t;

  tree_walk_cache_estep_baz(const unsigned b,
                            const unsigned s,
                            const rate_gtor_nscodon_base& _rg)
    : base_t(b,s), rg(_rg), int_node(0), out_node(0), bc() { }


  virtual
  void
  increment_site_stats(const tree_lhood_info& ti,
                       const site_data_fastlup& sdf,
                       const int site_id,
                       const site_data_fastlup_core::index_type* seq_state){
    local_increment_site_stats(ti,sdf,*this,site_id,seq_state,this->bc);
  }

  const rate_gtor_nscodon_base& rg;
  bi_tree_node* int_node;
  bi_tree_node* out_node;
  baz_counts bc;
};





static
void
local_increment_site_stats(const tree_lhood_info& ti,
                           const site_data_fastlup& sdf,
                           const tree_walk_cache_estep_baz& twe,
                           const int site_id,
                           const site_data_fastlup_core::index_type* seq_state,
                           baz_counts& bc){

  const smlfloat pnorm(static_cast<smlfloat>(sdf.data[site_id].unassigned_cat_count)*(1./twe.site_prob));

#if 0
  // do the root first
  for(int i(0);i<ti.n_states;++i){
    edat.root[i] +=
      pnorm*ti.tprob.root_prob[i]*
      twe.up_branch_prob.val[ti.tree.root->child1->node_index-1][i]*
      twe.up_branch_prob.val[ti.tree.root->child2->node_index-1][i];
  }
#endif

  // do the interior node (assuming a tree org topology) -- report for
  // each site

  if( static_cast<unsigned>(ti.n_states) != NSC4::SIZE ) { abort(); }

  prob_t* int_node_p(new prob_t[NSC4::SIZE]);
  prob_t* int_node_p_cdn(new prob_t[NSCODON::SIZE]);


  for(int i(0);i<ti.n_states;++i){
    int_node_p[i] =
      pnorm*
      twe.down_branch_prob.val[twe.int_node->node_index-1][i]*
      twe.up_branch_prob.val[twe.int_node->child1->node_index-1][i]*
      twe.up_branch_prob.val[twe.int_node->child2->node_index-1][i];
  }

  smlfloat sum = array_sum(int_node_p,NSC4::SIZE);
  pdistro_norm(int_node_p,int_node_p+NSC4::SIZE);

  // custom crap which assumes c4pre site
  nsc4_pdf_2_nscodon_pdf(int_node_p,int_node_p_cdn);
  for(unsigned i(0);i<NSCODON::SIZE;++i) int_node_p_cdn[i] *= sum;

  // deal with leaf branches separately:
  const int leaf1_id(twe.int_node->child1->leaf_index);
  const int leaf2_id(twe.int_node->child2->leaf_index);
  const int out_id(twe.out_node->leaf_index);
  if(leaf1_id >= 0 && leaf2_id >= 0){

    const unsigned leaf1_state(seq_state[leaf1_id]);
    const unsigned leaf2_state(seq_state[leaf2_id]);
    const unsigned out_state(seq_state[out_id]);

    const NSCODON::index_t lc1(NSC4::decode_nscodon(leaf1_state));
    const NSCODON::index_t lc2(NSC4::decode_nscodon(leaf2_state));
    const NSCODON::index_t lout(NSC4::decode_nscodon(out_state));

    // check lc1/lc2 separation conditions:
    NUC::index_t lc1n[CODON::BASE_SIZE];
    NUC::index_t lc2n[CODON::BASE_SIZE];
    NSCODON::decode(lc1n,static_cast<NSCODON::index_t>(lc1));
    NSCODON::decode(lc2n,static_cast<NSCODON::index_t>(lc2));


    int ndiff(0);
    for(unsigned n(0);n<CODON::BASE_SIZE;++n) if(lc1n[n] != lc2n[n]) ndiff++;


    const NSAA::index_t lc1aa(codon_trans_known(static_cast<NSCODON::index_t>(lc1)));
    const NSAA::index_t lc2aa(codon_trans_known(static_cast<NSCODON::index_t>(lc2)));

    if(lc1aa == lc2aa) {
      if(ndiff==1){
        bc.one_syn_on_leaf2 += int_node_p_cdn[lc1];
        bc.one_syn_on_leaf1 += int_node_p_cdn[lc2];
      }
    } else {
      if (ndiff==1) {
        bc.one_on_leaf2 += int_node_p_cdn[lc1];
        bc.one_on_leaf1 += int_node_p_cdn[lc2];

        bc.Kn_single_sub[lc1aa+NSAA::SIZE*lc2aa] += int_node_p_cdn[lc2];
        bc.Kn_single_sub[lc2aa+NSAA::SIZE*lc1aa] += int_node_p_cdn[lc1];

      // valid case to observe baz effect:
      } else if (ndiff==2) {
        // check that paths between codons are nonsyn:
        NSCODON::index_t i1,i2;
        if(lc1n[0] != lc2n[0]) {
          i1 = NSCODON::encode(lc1n[0],lc2n[1],lc2n[2]);
          i2 = NSCODON::encode(lc2n[0],lc1n[1],lc1n[2]);
        } else {
          i1 = NSCODON::encode(lc2n[0],lc1n[1],lc2n[2]);
          i2 = NSCODON::encode(lc1n[0],lc2n[1],lc1n[2]);
        }
        const NSAA::index_t i1aa(codon_trans_known(i1));
        const NSAA::index_t i2aa(codon_trans_known(i2));

        const bool is_2_nonsyn_away((i1aa != lc1aa && i1aa != lc2aa) &&
                                    (i2aa != lc1aa && i2aa != lc2aa));

        bool is_sneaky_path(false);

        if(is_2_nonsyn_away){
          // check against syn->non-syn->syn substitutions
          //
          // check among all synonymous mutations that are also two mutations away
          // thus, initial syn mut must be at one of the two differing nucleotides
          //
          for(unsigned nbase(0);nbase<CODON::BASE_SIZE;++nbase){
            if(lc1n[nbase] != lc2n[nbase]) {
              for(unsigned n(0);n<NUC::SIZE;++n){
                NUC::index_t nuc(static_cast<NUC::index_t>(n));
                if(nuc != lc1n[nbase] && nuc != lc2n[nbase]){
                  NUC::index_t ntmp[CODON::BASE_SIZE];
                  for(unsigned i(0);i<CODON::BASE_SIZE;++i) ntmp[i] = lc1n[i];
                  ntmp[nbase] = static_cast<NUC::index_t>(n);
                  const NSCODON::index_t testcodon(NSCODON::encode(ntmp));
                  const NSAA::index_t testaa(codon_trans_known(testcodon));

                  if(testaa == lc1aa){
                    NSCODON::index_t testi1,testi2;
                    if(ntmp[0] != lc2n[0]) {
                      testi1 = NSCODON::encode(ntmp[0],lc2n[1],lc2n[2]);
                      testi2 = NSCODON::encode(lc2n[0],ntmp[1],ntmp[2]);
                    } else {
                      testi1 = NSCODON::encode(lc2n[0],ntmp[1],lc2n[2]);
                      testi2 = NSCODON::encode(ntmp[0],lc2n[1],ntmp[2]);
                    }

                    const NSAA::index_t testi1aa(codon_trans_known(testi1));
                    const NSAA::index_t testi2aa(codon_trans_known(testi2));

                    if( ! ((testi1aa != lc1aa && testi1aa != lc2aa) &&
                           (testi2aa != lc1aa && testi2aa != lc2aa))){
                      is_sneaky_path=true;
                    }
                  }
                }
              }
            }
          }
        }


        if(is_2_nonsyn_away && (!is_sneaky_path)){

          // count 2 changes on one branch cases
          bc.bp[BP_BASIC].increment(int_node_p_cdn,lc1,lc2,i1,i2);

          bc.two_on_leaf_XX[lc1aa+NSAA::SIZE*lc2aa] += (int_node_p_cdn[lc1]+int_node_p_cdn[lc2])/2.;
          bc.two_on_leaf_XX[lc2aa+NSAA::SIZE*lc1aa] += (int_node_p_cdn[lc1]+int_node_p_cdn[lc2])/2.;

          bc.one_on_each_XX[lc1aa+NSAA::SIZE*lc2aa] += int_node_p_cdn[i1]+int_node_p_cdn[i2];
          bc.one_on_each_XX[lc2aa+NSAA::SIZE*lc1aa] += int_node_p_cdn[i1]+int_node_p_cdn[i2];


          if((lc1aa==NSAA::A && lc2aa==NSAA::R) ||
             (lc1aa==NSAA::R && lc2aa==NSAA::A)){
            bc.bp[BP_AR].increment(int_node_p_cdn,lc1,lc2,i1,i2);
          }

          if((lc1aa==NSAA::A && lc2aa==NSAA::I) ||
             (lc1aa==NSAA::I && lc2aa==NSAA::A)){
            bc.bp[BP_AI].increment(int_node_p_cdn,lc1,lc2,i1,i2);
           }

          if((lc1aa==NSAA::A && lc2aa==NSAA::N) ||
             (lc1aa==NSAA::N && lc2aa==NSAA::A)){
            bc.bp[BP_AN].increment(int_node_p_cdn,lc1,lc2,i1,i2);
          }

          if((is_hphobic(lc1aa) && is_polar(lc2aa)) ||
             (is_polar(lc1aa) && is_hphobic(lc2aa))){
            bc.bp[BP_HP].increment(int_node_p_cdn,lc1,lc2,i1,i2);
          }

          if(is_hphobic(lc1aa) && is_hphobic(lc2aa)){
            bc.bp[BP_HH].increment(int_node_p_cdn,lc1,lc2,i1,i2);
          }

          if(is_polar(lc1aa) && is_polar(lc2aa)){
            bc.bp[BP_PP].increment(int_node_p_cdn,lc1,lc2,i1,i2);
          }

          { // leaf 1
            NSAA::index_t highpath_aa(i1aa);
            NSAA::index_t lowpath_aa(i2aa);

            die("no category handlers in baz hack code");
            smlfloat i1score = std::sqrt(twe.rg.get_category_nsaa_param(lc2aa,i1aa,0,0)*
                                         twe.rg.get_category_nsaa_param(i1aa,lc1aa,0,0));

            smlfloat i2score = std::sqrt(twe.rg.get_category_nsaa_param(lc2aa,i2aa,0,0)*
                                         twe.rg.get_category_nsaa_param(i2aa,lc1aa,0,0));

            if( i1score < i2score ) {
              highpath_aa = i2aa;
              lowpath_aa = i1aa;
            }

            die("unhandled site/group categories");
            smlfloat val1 = twe.rg.get_category_nsaa_param(lc2aa,highpath_aa,0,0)*(int_node_p_cdn[lc2]);
            smlfloat val2 = twe.rg.get_category_nsaa_param(highpath_aa,lc1aa,0,0)*(int_node_p_cdn[lc2]);

            bc.Kn_leaf1_highpath_1 += val1;
            bc.Kn_leaf1_highpath_2 += val2;

            die("unhandled site/group categories");
            bc.Kn_leaf1_lowpath_1 += twe.rg.get_category_nsaa_param(lc2aa,lowpath_aa,0,0)*(int_node_p_cdn[lc2]);
            bc.Kn_leaf1_lowpath_2 += twe.rg.get_category_nsaa_param(lowpath_aa,lc1aa,0,0)*(int_node_p_cdn[lc2]);

            if(val1>val2){
              bc.Kn_leaf1_highpath_top += val1;
              bc.Kn_leaf1_highpath_bot += val2;
              bc.Kn_leaf1_highpath_top_frac += int_node_p_cdn[lc2];

              bc.Kn_highpath_top_sub[highpath_aa+NSAA::SIZE*lc2aa] += int_node_p_cdn[lc2];
              bc.Kn_highpath_bot_sub[lc1aa+NSAA::SIZE*highpath_aa] += int_node_p_cdn[lc2];
            } else {
              bc.Kn_leaf1_highpath_top += val2;
              bc.Kn_leaf1_highpath_bot += val1;

              bc.Kn_highpath_bot_sub[highpath_aa+NSAA::SIZE*lc2aa] += int_node_p_cdn[lc2];
              bc.Kn_highpath_top_sub[lc1aa+NSAA::SIZE*highpath_aa] += int_node_p_cdn[lc2];
            }

            bc.Kn_highpath_both_sub[lc1aa+NSAA::SIZE*lc2aa] += int_node_p_cdn[lc2];

            if(lc2aa==NSAA::A && lc1aa==NSAA::R){
              const NSCODON::index_t lc1i(static_cast<NSCODON::index_t>(lc1));
              if(lc1i==NSCODON::AGA || lc1i==NSCODON::AGG){
                bc.Kn_highpath_both_sub_ARlow += int_node_p_cdn[lc2];
              } else {
                bc.Kn_highpath_both_sub_ARhigh += int_node_p_cdn[lc2];
              }
            }

            bc.Kn_leaf1_highpath_diff += std::fabs(val1-val2);
          }


          { // leaf 2
            NSAA::index_t highpath_aa(i1aa);
            NSAA::index_t lowpath_aa(i2aa);

            smlfloat i1score = std::sqrt(twe.rg.get_category_nsaa_param(lc1aa,i1aa,0,0)*
                                         twe.rg.get_category_nsaa_param(i1aa,lc2aa,0,0));

            smlfloat i2score = std::sqrt(twe.rg.get_category_nsaa_param(lc1aa,i2aa,0,0)*
                                         twe.rg.get_category_nsaa_param(i2aa,lc2aa,0,0));

            if( i1score < i2score ) {
              highpath_aa = i2aa;
              lowpath_aa = i1aa;
            }

            smlfloat val1 = twe.rg.get_category_nsaa_param(lc1aa,highpath_aa,0,0)*(int_node_p_cdn[lc1]);
            smlfloat val2 = twe.rg.get_category_nsaa_param(highpath_aa,lc2aa,0,0)*(int_node_p_cdn[lc1]);

            bc.Kn_leaf2_highpath_1 += val1;
            bc.Kn_leaf2_highpath_2 += val2;

            bc.Kn_leaf2_lowpath_1 += twe.rg.get_category_nsaa_param(lc1aa,lowpath_aa,0,0)*(int_node_p_cdn[lc1]);
            bc.Kn_leaf2_lowpath_2 += twe.rg.get_category_nsaa_param(lowpath_aa,lc2aa,0,0)*(int_node_p_cdn[lc1]);

            if(val1>val2){
              bc.Kn_leaf2_highpath_top += val1;
              bc.Kn_leaf2_highpath_bot += val2;
              bc.Kn_leaf2_highpath_top_frac += int_node_p_cdn[lc1];
            } else {
              bc.Kn_leaf2_highpath_top += val2;
              bc.Kn_leaf2_highpath_bot += val1;

              bc.Kn_highpath_bot_sub[highpath_aa+NSAA::SIZE*lc1aa] += int_node_p_cdn[lc1];
              bc.Kn_highpath_top_sub[lc2aa+NSAA::SIZE*highpath_aa] += int_node_p_cdn[lc1];
            }

            bc.Kn_highpath_both_sub[lc2aa+NSAA::SIZE*lc1aa] += int_node_p_cdn[lc1];

            if(lc1aa==NSAA::A && lc2aa==NSAA::R){
              const NSCODON::index_t lc2i(static_cast<NSCODON::index_t>(lc2));
              if(lc2i==NSCODON::AGA || lc2i==NSCODON::AGG){
                bc.Kn_highpath_both_sub_ARlow += int_node_p_cdn[lc1];
              } else {
                bc.Kn_highpath_both_sub_ARhigh += int_node_p_cdn[lc1];
              }
            }

            bc.Kn_leaf2_highpath_diff += std::fabs(val1-val2);
          }


          bc.Kn_leaf1_1 += (twe.rg.get_category_nsaa_param(lc2aa,i1aa,0,0)+
                            twe.rg.get_category_nsaa_param(lc2aa,i2aa,0,0))*0.5*(int_node_p_cdn[lc2]);
          bc.Kn_leaf1_2 += (twe.rg.get_category_nsaa_param(i1aa,lc1aa,0,0)+
                            twe.rg.get_category_nsaa_param(i2aa,lc1aa,0,0))*0.5*(int_node_p_cdn[lc2]);
          bc.Kn_leaf2_1 += (twe.rg.get_category_nsaa_param(lc1aa,i1aa,0,0)+
                            twe.rg.get_category_nsaa_param(lc1aa,i2aa,0,0))*0.5*(int_node_p_cdn[lc1]);
          bc.Kn_leaf2_2 += (twe.rg.get_category_nsaa_param(i1aa,lc2aa,0,0)+
                            twe.rg.get_category_nsaa_param(i2aa,lc2aa,0,0))*0.5*(int_node_p_cdn[lc1]);

          bc.Kn_one_on_1 += (twe.rg.get_category_nsaa_param(i1aa,lc1aa,0,0)*(int_node_p_cdn[i1])+
                             twe.rg.get_category_nsaa_param(i2aa,lc1aa,0,0)*(int_node_p_cdn[i2]));
          bc.Kn_one_on_2 += (twe.rg.get_category_nsaa_param(i1aa,lc2aa,0,0)*(int_node_p_cdn[i1])+
                             twe.rg.get_category_nsaa_param(i2aa,lc2aa,0,0)*(int_node_p_cdn[i2]));


          if        (lout == lc1){
            bc.two_on_leaf2_parsi += 1;
          } else if (lout == lc2){
            bc.two_on_leaf1_parsi += 1;
          } else if (lout == i1 || lout == i2){
            bc.one_on_each_parsi += 1;
          }

          if(lc1n[1] == lc2n[1]){
            bc.bp[BP_13].increment(int_node_p_cdn,lc1,lc2,i1,i2);
          }
        }
      }
    }

  } else {
    die("unexpected option in baz");
  }

  delete [] int_node_p;
  delete [] int_node_p_cdn;
}




static
void
report_nsaa_grouping(const smlfloat* nsaa_selection,
                     const smlfloat* nsaa_flux,
                     const char* label,
                     const unsigned* reduction_map,
                     const char* syms,
                     const unsigned N,
                     std::ostream& os){

  smlfloat* selection_avg = new smlfloat[N*N];
  matrix_state_reduction(selection_avg,N,nsaa_selection,NSAA::SIZE,reduction_map,true);

  os << "MATRIX :: AVERAGE BY " << label << ":\n";
  matrix_report(selection_avg,N,syms,os);

  delete [] selection_avg;

  smlfloat* flux_sum = new smlfloat[N*N];
  matrix_state_reduction(flux_sum,N,nsaa_flux,NSAA::SIZE,reduction_map);

  os << "MATRIX :: SUM BY " << label << ":\n";
  matrix_report(flux_sum,N,syms,os);

  delete [] flux_sum;
}




static
void
dump_groupings(const char* label,
               prob_t* nsaa_matrix,
               std::ostream& os){

  { // report sub counts from aa-groupings
    using namespace NSAA_MAPS;

    os << label << " reduction groupings:\n";

    report_nsaa_grouping(nsaa_matrix,nsaa_matrix,"LEHNINGER GROUPS",
                         LEHNINGER::get_reduction_map(),LEHNINGER::syms,LEHNINGER::SIZE,os);

    report_nsaa_grouping(nsaa_matrix,nsaa_matrix,"LANGE GROUPS",
                         LANGE::get_reduction_map(),LANGE::syms,LANGE::SIZE,os);

    report_nsaa_grouping(nsaa_matrix,nsaa_matrix,"HP",
                         HP::get_reduction_map(),HP::syms,HP::SIZE,os);

    report_nsaa_grouping(nsaa_matrix,nsaa_matrix,"BETA BRANCHING",
                         BETA_BRANCH::get_reduction_map(),BETA_BRANCH::syms,BETA_BRANCH::SIZE,os);

    report_nsaa_grouping(nsaa_matrix,nsaa_matrix,"HELIX FORMING",
                         HELIX_FORM::get_reduction_map(),HELIX_FORM::syms,HELIX_FORM::SIZE,os);

    report_nsaa_grouping(nsaa_matrix,nsaa_matrix,"AROMATIC",
                         AROMATIC::get_reduction_map(),AROMATIC::syms,AROMATIC::SIZE,os);

    report_nsaa_grouping(nsaa_matrix,nsaa_matrix,"CYS",
                         CYS::get_reduction_map(),CYS::syms,CYS::SIZE,os);
    report_nsaa_grouping(nsaa_matrix,nsaa_matrix,"GLY",
                         GLY::get_reduction_map(),GLY::syms,GLY::SIZE,os);
    report_nsaa_grouping(nsaa_matrix,nsaa_matrix,"PRO",
                         PRO::get_reduction_map(),PRO::syms,PRO::SIZE,os);

    os << "\n\n";
  }
}






static
void
dump_doublemut_aa_paths(const rate_gtor_nscodon_base& rg,
                        std::ostream& os){


  for(unsigned c1(0);c1<NSCODON::SIZE;++c1){
    const NSCODON::index_t c1i(static_cast<NSCODON::index_t>(c1));
    NUC::index_t c1n1,c1n2,c1n3;
    NSCODON::decode(c1n1,c1n2,c1n3,c1i);
    const NSAA::index_t c1aa(codon_trans_known(c1i));

    for(unsigned c2(c1+1);c2<NSCODON::SIZE;++c2){
      const NSCODON::index_t c2i(static_cast<NSCODON::index_t>(c2));
      NUC::index_t c2n1,c2n2,c2n3;
      NSCODON::decode(c2n1,c2n2,c2n3,c2i);
      const NSAA::index_t c2aa(codon_trans_known(c2i));

      // check lc1/lc2 separation conditions:
      int ndiff(0);
      if(c1n1 != c2n1) ndiff++;
      if(c1n2 != c2n2) ndiff++;
      if(c1n3 != c2n3) ndiff++;

      if(c1aa != c2aa && ndiff==2) {

        // check that paths between codons are nonsyn:
        NSCODON::index_t i1,i2;
        if(c1n1 != c2n1) {
          i1 = NSCODON::encode(c1n1,c2n2,c2n3);
          i2 = NSCODON::encode(c2n1,c1n2,c1n3);
        } else {
          i1 = NSCODON::encode(c2n1,c1n2,c2n3);
          i2 = NSCODON::encode(c1n1,c2n2,c1n3);
        }
        const NSAA::index_t i1aa(codon_trans_known(static_cast<NSCODON::index_t>(i1)));
        const NSAA::index_t i2aa(codon_trans_known(static_cast<NSCODON::index_t>(i2)));

        if((i1aa != c1aa && i1aa != c2aa) &&
           (i2aa != c1aa && i2aa != c2aa)){

          { // leaf 1
            NSAA::index_t highpath_aa(i1aa);
            NSAA::index_t lowpath_aa(i2aa);

            smlfloat i1score = std::sqrt(rg.get_category_nsaa_param(c2aa,i1aa,0,0)*
                                         rg.get_category_nsaa_param(i1aa,c1aa,0,0));

            smlfloat i2score = std::sqrt(rg.get_category_nsaa_param(c2aa,i2aa,0,0)*
                                         rg.get_category_nsaa_param(i2aa,c1aa,0,0));

            if( i1score < i2score ) {
              highpath_aa = i2aa;
              lowpath_aa = i1aa;
            }

            os << "double_paths: leaf1 high: "
               << NSAA::syms[c2aa] << "->" << NSAA::syms[c1aa] << " "
               << NSAA::syms[c2aa] << "->" << NSAA::syms[highpath_aa] << "->" << NSAA::syms[c1aa]
               << " " << rg.get_category_nsaa_param(c2aa,highpath_aa,0,0)
               << " " << rg.get_category_nsaa_param(highpath_aa,c1aa,0,0) << "\n";

            os << "double_paths: leaf1 low_: "
               << NSAA::syms[c2aa] << "->" << NSAA::syms[c1aa] << " "
               << NSAA::syms[c2aa] << "->" << NSAA::syms[lowpath_aa] << "->" << NSAA::syms[c1aa]
               << " " << rg.get_category_nsaa_param(c2aa,lowpath_aa,0,0)
               << " " << rg.get_category_nsaa_param(lowpath_aa,c1aa,0,0) << "\n";
          }

          { // leaf 2
            NSAA::index_t highpath_aa(i1aa);
            NSAA::index_t lowpath_aa(i2aa);

            smlfloat i1score = std::sqrt(rg.get_category_nsaa_param(c1aa,i1aa,0,0)*
                                         rg.get_category_nsaa_param(i1aa,c2aa,0,0));

            smlfloat i2score = std::sqrt(rg.get_category_nsaa_param(c1aa,i2aa,0,0)*
                                         rg.get_category_nsaa_param(i2aa,c2aa,0,0));

            if( i1score < i2score ) {
              highpath_aa = i2aa;
              lowpath_aa = i1aa;
            }

            os << "double_paths: leaf2 high: "
               << NSAA::syms[c1aa] << "->" << NSAA::syms[c2aa] << " "
               << NSAA::syms[c1aa] << "->" << NSAA::syms[highpath_aa] << "->" << NSAA::syms[c2aa]
               << " " << rg.get_category_nsaa_param(c1aa,highpath_aa,0,0)
               << " " << rg.get_category_nsaa_param(highpath_aa,c2aa,0,0) << "\n";

            os << "double_paths: leaf2 low_: "
               << NSAA::syms[c1aa] << "->" << NSAA::syms[c2aa] << " "
               << NSAA::syms[c1aa] << "->" << NSAA::syms[lowpath_aa] << "->" << NSAA::syms[c2aa]
               << " " << rg.get_category_nsaa_param(c1aa,lowpath_aa,0,0)
               << " " << rg.get_category_nsaa_param(lowpath_aa,c2aa,0,0) << "\n";
          }
        }
      }
    }
  }

  os << "\n\n";
}
#endif // INCLUDE_BAZ_HACK




//entry point:
#ifndef INCLUDE_BAZ_HACK
void
lhood_site_cached_estep_baz(const rate_gtor_nscodon_base&,
                            const tree_lhood_info&,
                            const site_data_fastlup&){

  pass_away("lhood_site_cached_estep_baz: unmaintained -- assumes c4-pre site & no cats");
}
#else
void
lhood_site_cached_estep_baz(const rate_gtor_nscodon_base& rg,
                            const tree_lhood_info& ti,
                            const site_data_fastlup& sdf){

  site_data_fastlup_core::index_type* seq_state(new site_data_fastlup_core::index_type[ti.tree.n_leaves()]);

  int org_id(0);
  int site_id(0);

  tree_walk_cache_estep_baz twe(ti.tree.n_branches(),ti.n_states,rg);

  // find interior node
  if(ti.tree.root->child1->child1 != 0) {
    twe.int_node = ti.tree.root->child1;
    twe.out_node = ti.tree.root->child2;
  } else {
    twe.int_node = ti.tree.root->child2;
    twe.out_node = ti.tree.root->child1;
  }

  std::vector<tree_lhood_cache> prior_branch_lhood;
  lhood_site_cached_estep(ti,sdf,prior_branch_lhood,tree_lhood_cache(),org_id,seq_state,site_id,twe);

  delete [] seq_state;

  std::ostream& baz_fp(log_os);

  const char* c1_label(ti.tree.branch_label(twe.int_node->child1->node_index-1).c_str());
  const char* c2_label(ti.tree.branch_label(twe.int_node->child2->node_index-1).c_str());


  baz_fp << "leaf1 1sy " << c1_label << " " << twe.bc.one_syn_on_leaf1 << "\n";
  baz_fp << "leaf2 1sy " << c2_label << " " << twe.bc.one_syn_on_leaf2 << "\n";
  baz_fp << "\n";

  baz_fp << "leaf1 1ns " << c1_label << " " << twe.bc.one_on_leaf1 << "\n";
  baz_fp << "leaf2 1ns " << c2_label << " " << twe.bc.one_on_leaf2 << "\n";
  baz_fp << "\n";

  twe.bc.bp[BP_BASIC].report("",c1_label,c2_label,"xx",baz_fp);
  twe.bc.bp[BP_AR].report("AR",c1_label,c2_label,"xx",baz_fp);
  twe.bc.bp[BP_AI].report("AI",c1_label,c2_label,"xx",baz_fp);
  twe.bc.bp[BP_AN].report("AN",c1_label,c2_label,"xx",baz_fp);
  twe.bc.bp[BP_HP].report("HP",c1_label,c2_label,"xx",baz_fp);
  twe.bc.bp[BP_HH].report("HH",c1_label,c2_label,"xx",baz_fp);
  twe.bc.bp[BP_PP].report("PP",c1_label,c2_label,"xx",baz_fp);
  twe.bc.bp[BP_13].report("13",c1_label,c2_label,"xx",baz_fp);

  baz_fp << "leaf1 2ns_parsi " << c1_label << " "
            << twe.bc.two_on_leaf1_parsi << "\n";
  baz_fp << "leaf2 2ns_parsi " << c2_label << " "
            << twe.bc.two_on_leaf2_parsi << "\n";
  baz_fp << "each_ 2ns_parsi " << "xx " << twe.bc.one_on_each_parsi << "\n";
  baz_fp << "\n";


  baz_fp << "leaf1 2ns Kn1 " << c1_label << " " << twe.bc.Kn_leaf1_1/twe.bc.bp[BP_BASIC].two_on_leaf1 << "\n";
  baz_fp << "leaf1 2ns Kn2 " << c1_label << " " << twe.bc.Kn_leaf1_2/twe.bc.bp[BP_BASIC].two_on_leaf1 << "\n";
  baz_fp << "leaf2 2ns Kn1 " << c2_label << " " << twe.bc.Kn_leaf2_1/twe.bc.bp[BP_BASIC].two_on_leaf2 << "\n";
  baz_fp << "leaf2 2ns Kn2 " << c2_label << " " << twe.bc.Kn_leaf2_2/twe.bc.bp[BP_BASIC].two_on_leaf2 << "\n";


  baz_fp << "leaf1 2ns Kn1 high " << c1_label << " " << twe.bc.Kn_leaf1_highpath_1/twe.bc.bp[BP_BASIC].two_on_leaf1 << "\n";
  baz_fp << "leaf1 2ns Kn2 high " << c1_label << " " << twe.bc.Kn_leaf1_highpath_2/twe.bc.bp[BP_BASIC].two_on_leaf1 << "\n";
  baz_fp << "leaf2 2ns Kn1 high " << c2_label << " " << twe.bc.Kn_leaf2_highpath_1/twe.bc.bp[BP_BASIC].two_on_leaf2 << "\n";
  baz_fp << "leaf2 2ns Kn2 high " << c2_label << " " << twe.bc.Kn_leaf2_highpath_2/twe.bc.bp[BP_BASIC].two_on_leaf2 << "\n";


  baz_fp << "leaf1 2ns Kn1 low_ " << c1_label << " " << twe.bc.Kn_leaf1_lowpath_1/twe.bc.bp[BP_BASIC].two_on_leaf1 << "\n";
  baz_fp << "leaf1 2ns Kn2 low_ " << c1_label << " " << twe.bc.Kn_leaf1_lowpath_2/twe.bc.bp[BP_BASIC].two_on_leaf1 << "\n";
  baz_fp << "leaf2 2ns Kn1 low_ " << c2_label << " " << twe.bc.Kn_leaf2_lowpath_1/twe.bc.bp[BP_BASIC].two_on_leaf2 << "\n";
  baz_fp << "leaf2 2ns Kn2 low_ " << c2_label << " " << twe.bc.Kn_leaf2_lowpath_2/twe.bc.bp[BP_BASIC].two_on_leaf2 << "\n";


  baz_fp << "leaf1 2ns Kn high top " << c1_label << " " << twe.bc.Kn_leaf1_highpath_top/twe.bc.bp[BP_BASIC].two_on_leaf1 << "\n";
  baz_fp << "leaf1 2ns Kn high bot " << c1_label << " " << twe.bc.Kn_leaf1_highpath_bot/twe.bc.bp[BP_BASIC].two_on_leaf1 << "\n";
  baz_fp << "leaf2 2ns Kn high top " << c2_label << " " << twe.bc.Kn_leaf2_highpath_top/twe.bc.bp[BP_BASIC].two_on_leaf2 << "\n";
  baz_fp << "leaf2 2ns Kn high bot " << c2_label << " " << twe.bc.Kn_leaf2_highpath_bot/twe.bc.bp[BP_BASIC].two_on_leaf2 << "\n";

  baz_fp << "leaf1 2ns Kn high diff " << c1_label << " " << twe.bc.Kn_leaf1_highpath_diff/twe.bc.bp[BP_BASIC].two_on_leaf1 << "\n";
  baz_fp << "leaf2 2ns Kn high diff " << c2_label << " " << twe.bc.Kn_leaf2_highpath_diff/twe.bc.bp[BP_BASIC].two_on_leaf2 << "\n";

  baz_fp << "leaf1 2ns Kn high top frac " << c1_label << " " << twe.bc.Kn_leaf1_highpath_top_frac/twe.bc.bp[BP_BASIC].two_on_leaf1 << "\n";
  baz_fp << "leaf2 2ns Kn high top frac " << c2_label << " " << twe.bc.Kn_leaf2_highpath_top_frac/twe.bc.bp[BP_BASIC].two_on_leaf2 << "\n";


  baz_fp << "each_ 2ns 1 " << c1_label << " " << twe.bc.Kn_one_on_1/twe.bc.bp[BP_BASIC].one_on_each << "\n";
  baz_fp << "each_ 2ns 2 " << c2_label << " " << twe.bc.Kn_one_on_2/twe.bc.bp[BP_BASIC].one_on_each << "\n";
  baz_fp << "\n";


  // aa dependent ratios:
  //
  smlfloat XX_ratio[NSAA::SIZE*NSAA::SIZE];
  for(unsigned i(0);i<NSAA::SIZE*NSAA::SIZE;++i){
    if(twe.bc.two_on_leaf_XX[i] > 0.){
      XX_ratio[i] = twe.bc.two_on_leaf_XX[i]/twe.bc.one_on_each_XX[i];
    } else {
      XX_ratio[i] = 0.;
    }
  }

  for(unsigned i(0);i<NSAA::SIZE;++i){
    for(unsigned j(0);j<NSAA::SIZE;++j){
      if(i==j) continue;
      if(XX_ratio[j+NSAA::SIZE*i]==0.) continue;
      const unsigned joint(j+NSAA::SIZE*i);
      baz_fp << "split_ratio: " << NSAA::syms[i] << "->" << NSAA::syms[j]
             << " " << twe.bc.two_on_leaf_XX[joint] << " " << twe.bc.one_on_each_XX[joint]
             << " " << XX_ratio[joint] << "\n";

      baz_fp << "split_ratio: " << NSAA::syms[j] << "->" << NSAA::syms[i]
             << " " << twe.bc.two_on_leaf_XX[joint] << " " << twe.bc.one_on_each_XX[joint]
             << " " << XX_ratio[joint] << "\n";
    }
  }


  // normalize subs
  const smlfloat singlenorm(twe.bc.one_on_leaf1+twe.bc.one_on_leaf2);
  const smlfloat doublenorm(twe.bc.bp[BP_BASIC].two_on_leaf1+twe.bc.bp[BP_BASIC].two_on_leaf2);

  prob_t Kn_highpath_both_scaled[NSAA::SIZE*NSAA::SIZE];

  for(unsigned i(0);i<NSAA::SIZE*NSAA::SIZE;++i){
    twe.bc.Kn_single_sub[i] /= singlenorm;

    twe.bc.Kn_highpath_top_sub[i] /= doublenorm;
    twe.bc.Kn_highpath_bot_sub[i] /= doublenorm;

    Kn_highpath_both_scaled[i] = twe.bc.Kn_highpath_both_sub[i]/doublenorm;
  }

  twe.bc.Kn_highpath_both_sub_ARhigh /= doublenorm;
  twe.bc.Kn_highpath_both_sub_ARlow /= doublenorm;


  prob_t Kn_highpath_top_sub_norm[NSAA::SIZE*NSAA::SIZE];
  prob_t Kn_highpath_bot_sub_norm[NSAA::SIZE*NSAA::SIZE];

  // convert to ratios against single sub:
  for(unsigned i(0);i<NSAA::SIZE*NSAA::SIZE;++i){
    const smlfloat norm(twe.bc.Kn_single_sub[i]);
    if(norm != 0.) {
      Kn_highpath_top_sub_norm[i] = twe.bc.Kn_highpath_top_sub[i]/norm;
      Kn_highpath_bot_sub_norm[i] = twe.bc.Kn_highpath_bot_sub[i]/norm;
      //      twe.bc.Kn_highpath_both_sub[i] /= norm;
    } else {
      Kn_highpath_top_sub_norm[i] = 0;
      Kn_highpath_bot_sub_norm[i] = 0;
      //      twe.bc.Kn_highpath_both_sub[i] = 0;
    }
  }

  // --
  baz_fp << "\n\n";
  dump_groupings("SINGLE SUB",twe.bc.Kn_single_sub,baz_fp);

  // report top sub ratios:
  baz_fp << "\n\n";
  report_sorted_nsaa_pair_values(Kn_highpath_top_sub_norm,baz_fp,"HIGHPATH_TOP_NORM");
  baz_fp << "\n\n";

  dump_groupings("HIGHPATH TOP",twe.bc.Kn_highpath_top_sub,baz_fp);

  report_sorted_nsaa_pair_values(Kn_highpath_bot_sub_norm,baz_fp,"HIGHPATH_BOT_NORM");
  baz_fp << "\n\n";

  dump_groupings("HIGHPATH BOT",twe.bc.Kn_highpath_bot_sub,baz_fp);

  report_sorted_nsaa_pair_values(Kn_highpath_both_scaled,baz_fp,"HIGHPATH_BOTH");
  baz_fp << "\n\n";

  baz_fp << "HIGHPATH_BOTH A->Rhigh=" << twe.bc.Kn_highpath_both_sub_ARhigh << "\n";
  baz_fp << "HIGHPATH_BOTH A->Rlow=" << twe.bc.Kn_highpath_both_sub_ARlow << "\n";
  baz_fp << "\n\n";

  dump_groupings("HIGHPATH BOTH",Kn_highpath_both_scaled,baz_fp);
  baz_fp << "\n\n";

  dump_groupings("HIGHPATH BOTH RAW",twe.bc.Kn_highpath_both_sub,baz_fp);
  baz_fp << "\n\n";


  dump_doublemut_aa_paths(rg,baz_fp);
}
#endif // INCLUDE_BAZ_HACK


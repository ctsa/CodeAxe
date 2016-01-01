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

// $Id: rate_gtor_nscodon_base_report.cc 1217 2008-05-22 19:50:12Z ctsa $

/// \file

#include "bioseq_report.h"
#include "cat_manager.h"
#include "nsaa_maps.h"
#include "rate_gtor_nscodon_base.h"
#include "rate_gtor_nscodon_base_util.h"
#include "simple_util.h"
#include "stationary_pdistro.h"
#include "util/bio/bioseq_util_io.h"
#include "util/bio/bioseq_util_flux_rate_conversions.h"
#include "util/math/array_util.h"
#include "util/math/matrix_util.h"
#include "util/math/matrix_util_io.h"

#include <cmath>

#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>



const unsigned pdistro_row_write_width(7);



static
void
pdistro_row_write(const prob_t* pdistro,
                  const unsigned N,
                  std::ostream& os,
                  const unsigned width = pdistro_row_write_width,
                  const unsigned prec = 3){

  os.unsetf(std::ios::fixed | std::ios::scientific);
  os.setf(std::ios::showpoint);
  for(unsigned i(0);i<N;++i){
    os << " "
       << std::right
       << std::setw(width-1)
       << std::setprecision(prec)
       << pdistro[i];
  }
  os << "\n";
}



struct smc2{
  smc2(unsigned _i,
       unsigned _j,
       irv_t<smlfloat> _val1,
       irv_t<smlfloat> _val2,
       bool is_doswap=true)
  : i(_i),j(_j),val1(_val1),val2(_val2) {
    if(is_doswap){
      if(val1<val2){
        std::swap(val1,val2);
        std::swap(i,j);
      }
    }
    ratio = val1/val2;
  }

  bool
  operator<(const smc2& right) const {
    return ratio < right.ratio;
  }

  unsigned i;
  unsigned j;
  irv_t<smlfloat> val1;
  irv_t<smlfloat> val2;
  irv_t<smlfloat> ratio;
};


struct smc3{
  smc3(unsigned _i,
       unsigned _j,
       irv_t<smlfloat> _val1,
       irv_t<smlfloat> _val2)
  : i(_i),j(_j),val1(_val1),val2(_val2) { valdiff = val1-val2; }

  bool
  operator<(const smc3& right) const {
    return valdiff < right.valdiff;
  }

  unsigned i;
  unsigned j;
  irv_t<smlfloat> val1;
  irv_t<smlfloat> val2;
  irv_t<smlfloat> valdiff;
};




static
void
report_sorted_nsaa_exchange_parameters(const irv_t<smlfloat>* nsaa_rates,
                                       std::ostream& os,
                                       const char* label){

  report_sorted_nsaa_pair_values(nsaa_rates,os,label);

  {
    std::vector<smc2> vsm;

    /// dump rate pairs -- sorted by abs(log(a/b))
    for(unsigned i(0);i<NSAA::SIZE;++i){
      for(unsigned j(i+1);j<NSAA::SIZE;++j){
        if(i == j) continue;
        irv_t<smlfloat> val1 = nsaa_rates[j+i*NSAA::SIZE];
        if(std::fabs(val1) < std::numeric_limits<smlfloat>::epsilon() ) continue;
        irv_t<smlfloat> val2 = nsaa_rates[i+j*NSAA::SIZE];
        vsm.push_back(smc2(i,j,val1,val2));
      }
    }
    std::sort(vsm.begin(),vsm.end());
    std::reverse(vsm.begin(),vsm.end());

    os << "RATIO SORTED " << label << " PARAMETERS:\n";
    os << "(from->to=param & from->to=param & ratio )\n";
    for(unsigned i(0);i<vsm.size();++i){
      os << NSAA::syms[vsm[i].i] << "->"
         << NSAA::syms[vsm[i].j] << "="
         << std::setw(5)
         << std::left
         << vsm[i].val1 << " & "
         << NSAA::syms[vsm[i].j] << "->"
         << NSAA::syms[vsm[i].i] << "="
         << std::setw(5)
         << std::left
         << vsm[i].val2 << " & "
         << std::setw(4)
         << std::left
         << vsm[i].ratio;

      { // add codon transition list:
        std::vector<nscodon_pair> vc;
        get_codon_paths(vc,
                        static_cast<NSAA::index_t>(vsm[i].i),
                        static_cast<NSAA::index_t>(vsm[i].j));

        os << " &   "
           << NSAA::syms[vsm[i].i] << "->"
           << NSAA::syms[vsm[i].j] << " codons: ";
        for(unsigned j(0);j<vc.size();++j){
          if(j) os << ",";
          os << NSCODON::print(vc[j].c1) << "->" << NSCODON::print(vc[j].c2);
        }
      }
      os << "\n";
    }
    os << "\n";
  }
}




static
void
report_paired_aa_pdistro(const irv_t<smlfloat> distro1[NSAA::SIZE],
                         const irv_t<smlfloat> distro2[NSAA::SIZE],
                         std::ostream& os,
                         const char* label,
                         const char* label1,
                         const char* label2){

  //// dump sorted mut/acc rates:
  std::vector<std::pair<irv_t<smlfloat>,unsigned> > vs1,vs2;

  for(unsigned i(0);i<NSAA::SIZE;++i){
    vs1.push_back(std::make_pair(distro1[i],i));
    vs2.push_back(std::make_pair(distro2[i],i));
  }

  std::sort(vs1.begin(),vs1.end());
  std::reverse(vs1.begin(),vs1.end());
  std::sort(vs2.begin(),vs2.end());
  std::reverse(vs2.begin(),vs2.end());

  os << "SORTED " << label << " " << label1 << ":\n";
  for(unsigned i(0);i<NSAA::SIZE;++i)
    os << NSAA::syms[vs1[i].second] << "->X: "
       << std::setprecision(7) << vs1[i].first << "\n";
  os << "\n";

  os << "SORTED " << label << " " << label2 << ":\n";
  for(unsigned i(0);i<NSAA::SIZE;++i)
    os << "X->" << NSAA::syms[vs2[i].second] << ": "
       << std::setprecision(7) << vs2[i].first << "\n";
  os << "\n";


  std::vector<smc2> vs_ratio;
  for(unsigned i(0);i<NSAA::SIZE;++i){
    vs_ratio.push_back(smc2(i,i,distro1[i],distro2[i],false));
  }
  std::sort(vs_ratio.begin(),vs_ratio.end());
  std::reverse(vs_ratio.begin(),vs_ratio.end());

  os << label << " (" << label1 << "/" << label2 << ") RATIO SORTED:\n";
  os << "(from->to=param & from->to=param & ratio )\n";
  for(unsigned i(0);i<vs_ratio.size();++i){
    os << NSAA::syms[vs_ratio[i].i] << "->X="
       << std::setw(5)
       << std::left
       << vs_ratio[i].val1 << " & X->"
       << NSAA::syms[vs_ratio[i].i] << "="
       << std::setw(5)
       << std::left
       << vs_ratio[i].val2 << " & "
       << std::setw(4)
       << std::left
       << vs_ratio[i].ratio;
    os << "\n";
  }
  os << "\n";


  std::vector<smc3> vs_diff;
  for(unsigned i(0);i<NSAA::SIZE;++i){
    vs_diff.push_back(smc3(i,i,distro1[i],distro2[i]));
  }
  std::sort(vs_diff.begin(),vs_diff.end());
  std::reverse(vs_diff.begin(),vs_diff.end());

  os << label << " (" << label1 << " - " << label2 << ") DIFF SORTED:\n";
  os << "(from->to=param & from->to=param & diff(1-2) )\n";
  for(unsigned i(0);i<vs_diff.size();++i){
    os << NSAA::syms[vs_diff[i].i] << "->X="
       << std::setw(5)
       << std::left
       << vs_diff[i].val1 << " & X->"
       << NSAA::syms[vs_diff[i].i] << "="
       << std::setw(5)
       << std::left
       << vs_diff[i].val2 << " & "
       << std::setw(4)
       << std::left
       << vs_diff[i].valdiff;
    os << "\n";
  }
  os << "\n";
}




static
void
report_nsaa_exchange_averages(const irv_t<smlfloat> nsaa_mat[NSAA::SIZE*NSAA::SIZE],
                              std::ostream& os,
                              const char* label){

  // find simple averages for "to" and "from" values in matrix
  irv_t<smlfloat> from_avg[NSAA::SIZE];
  unsigned from_count[NSAA::SIZE];
  irv_t<smlfloat> to_avg[NSAA::SIZE];
  unsigned to_count[NSAA::SIZE];
  for(unsigned i(0);i<NSAA::SIZE;++i){
    from_avg[i] = 0.;
    from_count[i] = 0;
    to_avg[i] = 0.;
    to_count[i] = 0;
  }
  for(unsigned i(0);i<NSAA::SIZE;++i){
    for(unsigned j(0);j<NSAA::SIZE;++j){
      if(i == j) continue;
      const irv_t<smlfloat> val = nsaa_mat[j+i*NSAA::SIZE];
      if(std::fabs(val) < std::numeric_limits<smlfloat>::epsilon() ) continue;
      from_avg[i] += val;
      from_count[i]++;
      to_avg[j] += val;
      to_count[j]++;
    }
  }
  for(unsigned i(0);i<NSAA::SIZE;++i){
    from_avg[i] /= static_cast<smlfloat>(from_count[i]);
    to_avg[i] /= static_cast<smlfloat>(to_count[i]);
  }

  report_paired_aa_pdistro(from_avg,to_avg,os,label,"FROM AVG","TO AVG");
}




static
void
report_nsaa_exchange_sums(const irv_t<smlfloat> nsaa_mat[NSAA::SIZE*NSAA::SIZE],
                          std::ostream& os,
                          const char* label){

  // find simple averages for to and from values in matrix
  irv_t<smlfloat> from_stat[NSAA::SIZE];
  irv_t<smlfloat> to_stat[NSAA::SIZE];
  for(unsigned i(0);i<NSAA::SIZE;++i){
    from_stat[i] = 0.;
    to_stat[i] = 0.;
  }
  for(unsigned i(0);i<NSAA::SIZE;++i){
    for(unsigned j(0);j<NSAA::SIZE;++j){
      if(i == j) continue;
      const irv_t<smlfloat> val = nsaa_mat[j+i*NSAA::SIZE];
      if(std::fabs(val) < std::numeric_limits<smlfloat>::epsilon() ) continue;
      from_stat[i] += val;
      to_stat[j] += val;
    }
  }

  report_paired_aa_pdistro(from_stat,to_stat,os,label,"FROM ALL","TO ALL");
}




static
void
report_nsaa_matrix_full(const irv_t<smlfloat>* nsaa_mat,
                        std::ostream& os,
                        const char* label){

  os << label << " MATRIX:\n";
  matrix_report(nsaa_mat,NSAA::SIZE,NSAA::syms,os);

  report_sorted_nsaa_exchange_parameters(nsaa_mat,os,label);
}


#if 0
static
void
report_nsaa_grouping(const irv_t<smlfloat>* nsaa_selection,
                     const irv_t<smlfloat>* nsaa_flux,
                     const char* label,
                     const unsigned* reduction_map,
                     const char* syms,
                     const unsigned N,
                     std::ostream& os){

  irv_t<smlfloat>* selection_avg(new irv_t<smlfloat>[N*N]);
  matrix_state_reduction(selection_avg,N,nsaa_selection,NSAA::SIZE,reduction_map,true);

  os << "AA SELECTION :: AVERAGE BY " << label << ":\n";
  matrix_report(selection_avg,N,syms,os);

  irv_t<smlfloat>* flux_sum(selection_avg);
  matrix_state_reduction(flux_sum,N,nsaa_flux,NSAA::SIZE,reduction_map);

  os << "AA FLUX :: SUM BY " << label << ":\n";
  matrix_report(flux_sum,N,syms,os);

  delete [] flux_sum;
}
#endif





// extremely sloppy amino acid mutability set dfs
static
void
aaset_dfs(const bool* is_connected,
          const unsigned* aa_list,
          std::vector<bool>& hit_list,
          const unsigned N,
          const unsigned x){

  hit_list[x] = true;

  for(unsigned i(0);i<N;++i){
    if(i==x) continue;
    if( is_connected[aa_list[i]+NSAA::SIZE*aa_list[x]] ){
      if(! hit_list[i]){
        aaset_dfs(is_connected,aa_list,hit_list,N,i);
      }
    }
  }
}



static
bool
is_aaset_connected(const bool* is_connected,
                   const unsigned* aa_list,
                   const unsigned N){

  std::vector<bool> hit_list(N,false);

  // run a dfs, see if anything isn't colored
  aaset_dfs(is_connected,aa_list,hit_list,N,0);

  for(unsigned i(0);i<N;++i) if(! hit_list[i]) return false;

  return true;
}



static
irv_t<smlfloat>
get_Nset_selavg(const irv_t<smlfloat>* neutral_nsaa_flux,
                const irv_t<smlfloat>* nsaa_flux,
                const unsigned* aa_list,
                const unsigned N){

  irv_t<smlfloat> sum(0.);
  unsigned count(0);
  for(unsigned i(0);i<N;++i){
    for(unsigned j(0);j<N;++j){
      if(i==j) continue;
      if(neutral_nsaa_flux[aa_list[i]+aa_list[j]*NSAA::SIZE] <= 0.) continue;
      sum += nsaa_flux[aa_list[i]+aa_list[j]*NSAA::SIZE]/
        neutral_nsaa_flux[aa_list[i]+aa_list[j]*NSAA::SIZE];
      count++;
    }
  }
  if(count) return sum/static_cast<smlfloat>(count);
  else      return irv_t<smlfloat>(0.);
}


struct aaset_selavg_score {

  irv_t<smlfloat>
  score(unsigned* aaid_set,
        const unsigned set_size) const {
    return get_Nset_selavg(neutral_nsaa_flux,nsaa_flux,aaid_set,set_size);
  }

  const irv_t<smlfloat>* neutral_nsaa_flux;
  const irv_t<smlfloat>* nsaa_flux;
};



static
irv_t<smlfloat>
get_Nset_dn_ds(const irv_t<smlfloat>* neutral_nsaa_flux,
               const irv_t<smlfloat>* nsaa_flux,
               const unsigned* aa_list,
               const unsigned N){

  smlfloat flux (0.);
  smlfloat nflux(0.);
  for(unsigned i(0);i<N;++i){
    for(unsigned j(0);j<N;++j){
      if(i==j) continue;
      flux += nsaa_flux[aa_list[i]+aa_list[j]*NSAA::SIZE];
      nflux += neutral_nsaa_flux[aa_list[i]+aa_list[j]*NSAA::SIZE];
    }
  }
  if(nflux>0) return flux/nflux;
  else        return irv_t<smlfloat>(0.);
}


struct aaset_dn_ds_score {

  irv_t<smlfloat>
  score(unsigned* aaid_set,
        const unsigned set_size) const {
    return get_Nset_dn_ds(neutral_nsaa_flux,nsaa_flux,aaid_set,set_size);
  }

  const irv_t<smlfloat>* neutral_nsaa_flux;
  const irv_t<smlfloat>* nsaa_flux;
};



template <typename scorefunc>
void
score_aaset(const scorefunc& sf,
            const bool* is_reachable,
            unsigned* aaid_set,
            const unsigned set_size,
            unsigned s,
            std::ostream& os){

  const unsigned end(s==0 ? static_cast<unsigned>(NSAA::SIZE) : aaid_set[s-1]);

  for(unsigned i(0);i<end;++i){
    aaid_set[s] = i;
    if((s+1)==set_size){
      if(! is_aaset_connected(is_reachable,aaid_set,set_size)) continue;
      os << "GKEY " << set_size << " ";
      os << std::setw(9);
      os << sf.score(aaid_set,set_size);
      for(unsigned j(0);j<set_size;++j){
        os << " " << NSAA::syms[aaid_set[j]] << " ";
      }
      os << "\n";
    }else {
      score_aaset(sf,is_reachable,aaid_set,set_size,s+1,os);
    }
  }
}


#ifdef REPORT_CLUSTERS
template <typename scorefunc>
void
score_aaset_start(const irv_t<smlfloat>* neutral_nsaa_flux,
                  const irv_t<smlfloat>* nsaa_flux,
                  const bool* is_reachable,
                  std::ostream& os){

  const unsigned max_set_size(11);
  scorefunc sf;
  sf.neutral_nsaa_flux=neutral_nsaa_flux;
  sf.nsaa_flux=nsaa_flux;

  for(unsigned set_size(2);set_size<max_set_size;++set_size){
    simple_array<unsigned> aaid_set(set_size);
    score_aaset(sf,is_reachable,aaid_set.ptr(),set_size,0,os);
  }
}



static
void
report_nsaa_sets(const irv_t<smlfloat>* neutral_nsaa_flux,
                 const irv_t<smlfloat>* nsaa_flux,
                 std::ostream& os){

  //  irv_t<smlfloat>* pair_dn_ds(new irv_t<smlfloat>[NSAA::SIZE*NSAA::SIZE]);

  bool is_reachable[NSAA::SIZE*NSAA::SIZE];

  make_reachable_aa_matrix(is_reachable);
  os << std::setprecision(6);

  // use dn/ds sets
  //  score_aaset_start<aaset_dn_ds_score>(neutral_nsaa_flux,nsaa_flux,is_reachable,os);

  // use selection strength average sets...
  score_aaset_start<aaset_selavg_score>(neutral_nsaa_flux,nsaa_flux,is_reachable,os);
}
#endif  // REPORT_CLUSTERS



#if 0
  // over all pairs:
  for(unsigned i(0);i<NSAA::SIZE;++i){
    for(unsigned j(i+1);j<NSAA::SIZE;++j){

#if 0
      if(neutral_nsaa_flux[j+i*NSAA::SIZE] <= 0.) {
        pair_dn_ds[j+i*NSAA::SIZE]=0.;
      } else {
        smlfloat flux(nsaa_flux[j+i*NSAA::SIZE]+nsaa_flux[i+j*NSAA::SIZE]);
        smlfloat nflux(neutral_nsaa_flux[j+i*NSAA::SIZE]+neutral_nsaa_flux[i+j*NSAA::SIZE]);
        pair_dn_ds[j+i*NSAA::SIZE] = flux/nflux;
      }
#endif
#if 0
      for(unsigned k(j+1);k<NSAA::SIZE;++k){

        unsigned nsfoo[] = {i,j,k};
        if(! is_aaset_connected(is_reachable,nsfoo,3)) continue;
        os << "GKEY " << NSAA::syms[i] << " " << NSAA::syms[j] << " " << NSAA::syms[k] << " ";
        os << get_Nset_dn_ds(neutral_nsaa_flux,nsaa_flux,3,nsfoo);
        os << "\n";
      }
#endif
    }
  }
#if 0
  os << "PAIR DN/DS\n";
  matrix_report(pair_dn_ds,NSAA::SIZE,NSAA::syms,os);

  report_nsaa_matrix_full(pair_dn_ds,os,"PAIR DN/DS");
#endif
#endif




static
void
report_nsaa_grouping_flux_average(const irv_t<smlfloat>* neutral_nsaa_flux,
                                  const irv_t<smlfloat>* nsaa_flux,
                                  const char* label,
                                  const unsigned* reduction_map,
                                  const char* syms,
                                  const unsigned N,
                                  std::ostream& os){

  irv_t<smlfloat>* neutral_reduced_flux(new irv_t<smlfloat>[N*N]);
  matrix_state_reduction(neutral_reduced_flux,N,neutral_nsaa_flux,NSAA::SIZE,reduction_map);
#if 0
  os << "NEUTRAL AA FLUX :: SUM BY " << label << ":\n";
  matrix_report(neutral_reduced_flux,N,syms,os);
#endif
  irv_t<smlfloat>* reduced_flux(new irv_t<smlfloat>[N*N]);
  matrix_state_reduction(reduced_flux,N,nsaa_flux,NSAA::SIZE,reduction_map);

  irv_t<smlfloat>* reduced_dn_ds(new irv_t<smlfloat>[N*N]);

  for(unsigned i(0);i<(N*N);++i){
    if(neutral_reduced_flux[i]<=0.) reduced_dn_ds[i] = 0.;
    else                            reduced_dn_ds[i] = reduced_flux[i]/neutral_reduced_flux[i];
  }

  static const unsigned rwid(9);
  static const unsigned rpre(3);

  os << "AA SELECTION :: NSAA FLUX AVERAGE BY " << label << ":\n";
  matrix_report(reduced_dn_ds,N,syms,os,rwid,rpre);

  os << "AA FLUX :: SUM BY " << label << ":\n";
  matrix_report(reduced_flux,N,syms,os,rwid,rpre);

  os << "AA NEUTRAL FLUX :: SUM BY " << label << ":\n";
  matrix_report(neutral_reduced_flux,N,syms,os,rwid,rpre);

  delete [] neutral_reduced_flux;
  delete [] reduced_flux;
  delete [] reduced_dn_ds;
}



static
void
print_section_divider(std::ostream& os){
  os << "######################################################################\n";
}



static
void
print_branch_cat_sel_matrix_header(const char* subsection_label,
                                   const unsigned branch_set_id,
                                   const std::string& branch_set_label,
                                   const unsigned cat_id,
                                   const std::string& cat_label,
                                   const unsigned sel_matrix_model_id,
                                   const std::string& sel_matrix_model_label,
                                   std::ostream& os){
  print_section_divider(os);
  os << "branch_cat_set_seq_cat_sel_matrix_model_summary: " << subsection_label << "\n\n";
  os << "branch_cat_set: " << branch_set_id << " " << branch_set_label << "\n";
  os << "global_seq_cat: " << cat_id << " " << cat_label << "\n";
  os << "sel_matrix_model_param_cat: " << sel_matrix_model_id << " " << sel_matrix_model_label << "\n\n";
}



struct sel_cat_data {
  sel_cat_data()
    : prior(0.), exp_dn_ds(0.), label() {}

  sel_cat_data(irv_t<prob_t> p,
               irv_t<smlfloat> e,
               const std::string& l)
    : prior(p), exp_dn_ds(e), label(l) {}

  bool operator<(const sel_cat_data& rhs) const { return (rhs.prior < rhs.prior);}

  irv_t<prob_t> prior;
  irv_t<smlfloat> exp_dn_ds;
  std::string label;
};



void
rate_gtor_nscodon_base::
report(std::ostream& os) const {

  // report from mutation model:
  base_t::report(os);

  {
    using namespace CODON_BIAS_MODEL;

    const index_t& cbm(_data.copt.codon_bias_model);
    if(cbm == SYNON_RATIO ||
       cbm == SYNON_NORM_RATIO){
      // dump codon bias parameters
      os << "CODON BIAS SELECTION:\n";
      os << std::setprecision(6) << std::fixed;
      for(unsigned i(0);i<NSCODON::SIZE;++i){
        const NSCODON::index_t nsi(static_cast<NSCODON::index_t>(i));
        os << NSCODON::print(nsi) << " : " << codon_bias_factor(nsi) << "\n";
      }
      os << "\n";
    }
  }

  irv_t<smlfloat> nsaa_selection[NSAA::SIZE*NSAA::SIZE];
  irv_t<smlfloat> neutral_nscodon_flux[NSCODON::SIZE*NSCODON::SIZE];
  irv_t<smlfloat> neutral_nsaa_flux[NSAA::SIZE*NSAA::SIZE];
  irv_t<smlfloat> nsaa_flux[NSAA::SIZE*NSAA::SIZE];
  irv_t<smlfloat> nsaa_rate[NSAA::SIZE*NSAA::SIZE];
  prob_t bg_pdistro_nsaa_cat[NSAA::SIZE];

  const unsigned nsmc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_MATRIX));
  const unsigned nssc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_STRENGTH,CAT_MIX_TYPE::SITE));
  const unsigned ngsc(cm().typed_cat_size(CAT_PARAM_TYPE::SEL_STRENGTH,CAT_MIX_TYPE::GROUP));
  const unsigned n_cats(cm().cat_size());

  simple_array<prob_t> cprob(n_cats);
  cm().cat_pdistro(cprob.ptr());

  /// calculate expected dN/dS for each category on each branch
  ///
  /// \todo make this report more efficient by grouping together branches with the same
  /// branch category pattern across all cats...
  ///
  const unsigned n_branches(cm().branch_size());

  for(unsigned b(0);b<n_branches;++b){

    const std::string& branch_label(cm().branch_label(b));

    irv_t<smlfloat> exp_dn_ds(0.);
    std::vector<irv_t<smlfloat> > smc_exp_dn_ds(nsmc,0.);
    std::vector<irv_t<smlfloat> > ssc_exp_dn_ds(nssc,0.);
    std::vector<irv_t<smlfloat> > gsc_exp_dn_ds(ngsc,0.);

    {
      std::vector<irv_t<smlfloat> > smc_exp_dn_ds_norm(nsmc,0.);
      std::vector<irv_t<smlfloat> > ssc_exp_dn_ds_norm(nssc,0.);
      std::vector<irv_t<smlfloat> > gsc_exp_dn_ds_norm(ngsc,0.);

      for(unsigned c(0);c<n_cats;++c){
        const unsigned branch_cat_set(cm().get_branch_cat_set(c,b));

        // get dn/ds(cat,branch):
        nscodon_flux_no_aa_selection(*this,c,branch_cat_set,neutral_nscodon_flux);
        nsaa_selection_matrix(nsaa_selection,c,branch_cat_set);
        const smlfloat cat_branch_dn_ds(neutral_nscodon_flux_to_dn_ds(neutral_nscodon_flux,nsaa_selection));

        const unsigned smc(cm().typed_cat_no_from_cat_no_y_branch_id(c,b,CAT_PARAM_TYPE::SEL_MATRIX));
        const unsigned ssc(cm().typed_cat_no_from_cat_no_y_branch_id(c,b,CAT_PARAM_TYPE::SEL_STRENGTH,CAT_MIX_TYPE::SITE));
        const unsigned gsc(cm().typed_cat_no_from_cat_no_y_branch_id(c,b,CAT_PARAM_TYPE::SEL_STRENGTH,CAT_MIX_TYPE::GROUP));

        const prob_t cp(cprob[c]);

        smc_exp_dn_ds[smc] += cat_branch_dn_ds*cp;
        smc_exp_dn_ds_norm[smc] += cp;
        ssc_exp_dn_ds[ssc] += cat_branch_dn_ds*cp;
        ssc_exp_dn_ds_norm[ssc] += cp;
        gsc_exp_dn_ds[gsc] += cat_branch_dn_ds*cp;
        gsc_exp_dn_ds_norm[gsc] += cp;
        exp_dn_ds += cat_branch_dn_ds*cp;
      }

      for(unsigned smc(0);smc<nsmc;++smc){ smc_exp_dn_ds[smc] /= smc_exp_dn_ds_norm[smc]; }
      for(unsigned ssc(0);ssc<nssc;++ssc){ ssc_exp_dn_ds[ssc] /= ssc_exp_dn_ds_norm[ssc]; }
      for(unsigned gsc(0);gsc<ngsc;++gsc){ gsc_exp_dn_ds[gsc] /= gsc_exp_dn_ds_norm[gsc]; }
    }

    print_section_divider(os);
    os << "Reporting selection for branch [" << branch_label << "]\n\n";

    os.unsetf(std::ios::fixed);
    os << std::setprecision(5);
    os << std::left;
    os << "branch_expected_dN/dS: " << exp_dn_ds << "\n\n";

    {
      // dN/dS for select matrix cats:
      //
      simple_array<prob_t> smc_prior(nsmc);
      cm().branch_typed_cat_pdistro(smc_prior.ptr(),b,CAT_PARAM_TYPE::SEL_MATRIX);

      os << "SITE SELECTION MATRIX CATEGORIES (" << nsmc << ") dN/dS:\n";

      for(unsigned smc(0);smc<nsmc;++smc){
        os << "[ " <<  cm().typed_cat_label(smc,CAT_PARAM_TYPE::SEL_MATRIX) << " ]"
           << " branch_expected_dN/dS= " << std::setw(7) << smc_exp_dn_ds[smc]
           << " branch_prior_prob= " << std::setw(7) << smc_prior[smc] << "\n";
      }
      os << "\n";
    }

    {
      // dN/dS for site selection cats:
      //
      simple_array<prob_t> ssc_prior(nssc);
      cm().branch_typed_cat_pdistro(ssc_prior.ptr(),b,CAT_PARAM_TYPE::SEL_STRENGTH,CAT_MIX_TYPE::SITE);

      std::vector<std::pair<irv_t<smlfloat>,sel_cat_data> > tmp_cats(nssc);
      for(unsigned ssc(0);ssc<nssc;++ssc){
        const std::string& cat_label(cm().typed_cat_label(ssc,CAT_PARAM_TYPE::SEL_STRENGTH,CAT_MIX_TYPE::SITE));
        tmp_cats[ssc] = std::make_pair(get_paramset_member(PARAMSET_SITE_SELECT_CATS,ssc),
                                       sel_cat_data(ssc_prior[ssc],ssc_exp_dn_ds[ssc],cat_label));
      }
      std::sort(tmp_cats.begin(),tmp_cats.end());

      os << "SITE SELECTION CATEGORIES (" << nssc << "):\n";
      for(unsigned ssc(0);ssc<nssc;++ssc){
        os << "[ " << tmp_cats[ssc].second.label << " ]"
           << " scale_factor= " << std::setw(7) << tmp_cats[ssc].first
           << " branch_expected_dN/dS= " << std::setw(7) << tmp_cats[ssc].second.exp_dn_ds
           << " branch_prior_prob= " << std::setw(7) << tmp_cats[ssc].second.prior << "\n";
      }
      os << "\n";
    }

    {
      // dN/dS for group selection cats
      //
      simple_array<prob_t> gsc_prior(ngsc);
      cm().branch_typed_cat_pdistro(gsc_prior.ptr(),b,CAT_PARAM_TYPE::SEL_STRENGTH,CAT_MIX_TYPE::GROUP);

      std::vector<std::pair<irv_t<smlfloat>,sel_cat_data> > tmp_cats(ngsc);
      for(unsigned i(0);i<ngsc;++i){
        const std::string& cat_label(cm().typed_cat_label(i,CAT_PARAM_TYPE::SEL_STRENGTH,CAT_MIX_TYPE::GROUP));
        tmp_cats[i] = std::make_pair(get_paramset_member(PARAMSET_GROUP_SELECT_CATS,i),
                                     sel_cat_data(gsc_prior[i],gsc_exp_dn_ds[i],cat_label));
      }
      std::sort(tmp_cats.begin(),tmp_cats.end());

      os << "GROUP SELECTION CATEGORIES (" << ngsc << "):\n";
      for(unsigned i(0);i<ngsc;++i){
        os << "[ " << tmp_cats[i].second.label << " ]"
           << " scale_factor= " << std::setw(7) << tmp_cats[i].first
           << " branch_expected_dN/dS= " << std::setw(7) << tmp_cats[i].second.exp_dn_ds
           << " branch_prior_prob= " << std::setw(7) << tmp_cats[i].second.prior << " \n";
      }
      os << "\n";
    }
  }



  // it's no longer simple to determine when select matrices can be
  // summarized without taking into account the category background
  // the selection matrix is in, so the code now rolls through and
  // reports on all cats to insure correctness
  //
  for(unsigned c(0);c<n_cats;++c){
    const rates_func_options_base bopt(c);

    const unsigned n_branch_cat_sets(cm().branch_cat_set_size(c));

    for(unsigned bcs(0);bcs<n_branch_cat_sets;++bcs){
      std::string branch_set_label;
      cm().branch_cat_set_label(c,bcs,branch_set_label);

      {
        nscodon_flux_no_aa_selection(*this,c,bcs,neutral_nscodon_flux);
        nsaa_selection_matrix(nsaa_selection,c,bcs);
        nscodon_pdf_2_nsaa_pdf(bg_pdistro_nscodon_cat(c),bg_pdistro_nsaa_cat);

        const std::string& cat_label(cm().cat_label(c));

        const unsigned smc(cm().typed_cat_no_from_cat_no_y_branch_cat_set(c,bcs,CAT_PARAM_TYPE::SEL_MATRIX));

        if(select_model(smc) == SELECT_MODEL::NONE) continue;

        const std::string& mcat_label(cm().typed_cat_label(smc,CAT_PARAM_TYPE::SEL_MATRIX));
        print_branch_cat_sel_matrix_header("AA_SELECTION",bcs,branch_set_label,c,cat_label,smc,mcat_label,os);

        if(select_model(smc) == SELECT_MODEL::SINGLE) {
          const smlfloat cat_branch_dn_ds(neutral_nscodon_flux_to_dn_ds(neutral_nscodon_flux,nsaa_selection));
          os << "cat_branch_set_dN/dS: " << cat_branch_dn_ds << "\n\n";
          continue;
        }

        report_nsaa_matrix_full(nsaa_selection,os,"AA SELECTION");
        report_nsaa_exchange_averages(nsaa_selection,os,"AA SELECTION");

        {  // average aa selection from the background codon distro and neutral rates
          irv_t<smlfloat> to_avg[NSAA::SIZE];
          irv_t<smlfloat> from_avg[NSAA::SIZE];
          neutral_nscodon_flux_to_nsaa_dn_ds(neutral_nscodon_flux,nsaa_selection,to_avg,from_avg);

          report_paired_aa_pdistro(from_avg,to_avg,os,"AA SELECTION/NSCODON AVG","FROM AVG","TO AVG");
        }

        print_branch_cat_sel_matrix_header("AA_FLUX",bcs,branch_set_label,c,cat_label,smc,mcat_label,os);

        neutral_nscodon_flux_to_nsaa_flux(neutral_nscodon_flux,nsaa_selection,nsaa_flux);
        report_nsaa_matrix_full(nsaa_flux,os,"AA FLUX");
        report_nsaa_exchange_sums(nsaa_flux,os,"AA FLUX");


        print_branch_cat_sel_matrix_header("AA_RATE",bcs,branch_set_label,c,cat_label,smc,mcat_label,os);

        std::copy(nsaa_flux,nsaa_flux+NSAA::SIZE*NSAA::SIZE,nsaa_rate);

        flux_to_rates_inplace(NSAA::SIZE,bg_pdistro_nsaa_cat,nsaa_rate);

        report_nsaa_matrix_full(nsaa_rate,os,"AA RATE");

        nscodon_flux_to_nsaa_flux(neutral_nscodon_flux,neutral_nsaa_flux);

        { // report exchanges from aa-groupings
          using namespace NSAA_MAPS;
          print_branch_cat_sel_matrix_header("AA_SUMMARY",bcs,branch_set_label,c,cat_label,smc,mcat_label,os);

          report_nsaa_grouping_flux_average(neutral_nsaa_flux,nsaa_flux,"LEHNINGER GROUPS",
                                            LEHNINGER::get_reduction_map(),LEHNINGER::syms,LEHNINGER::SIZE,os);

          report_nsaa_grouping_flux_average(neutral_nsaa_flux,nsaa_flux,"LANGE GROUPS",
                                            LANGE::get_reduction_map(),LANGE::syms,LANGE::SIZE,os);

          report_nsaa_grouping_flux_average(neutral_nsaa_flux,nsaa_flux,"HP",
                                            HP::get_reduction_map(),HP::syms,HP::SIZE,os);

          report_nsaa_grouping_flux_average(neutral_nsaa_flux,nsaa_flux,"HPO",
                                            HPO::get_reduction_map(),HPO::syms,HPO::SIZE,os);

          report_nsaa_grouping_flux_average(neutral_nsaa_flux,nsaa_flux,"BETA BRANCHING",
                                            BETA_BRANCH::get_reduction_map(),BETA_BRANCH::syms,BETA_BRANCH::SIZE,os);

          report_nsaa_grouping_flux_average(neutral_nsaa_flux,nsaa_flux,"HELIX FORMING",
                                            HELIX_FORM::get_reduction_map(),HELIX_FORM::syms,HELIX_FORM::SIZE,os);

          report_nsaa_grouping_flux_average(neutral_nsaa_flux,nsaa_flux,"AROMATIC",
                                            AROMATIC::get_reduction_map(),AROMATIC::syms,AROMATIC::SIZE,os);

          report_nsaa_grouping_flux_average(neutral_nsaa_flux,nsaa_flux,"CYS",
                                            CYS::get_reduction_map(),CYS::syms,CYS::SIZE,os);

          report_nsaa_grouping_flux_average(neutral_nsaa_flux,nsaa_flux,"GLY",
                                            GLY::get_reduction_map(),GLY::syms,GLY::SIZE,os);

          report_nsaa_grouping_flux_average(neutral_nsaa_flux,nsaa_flux,"PRO",
                                            PRO::get_reduction_map(),PRO::syms,PRO::SIZE,os);
        }

#ifdef REPORT_CLUSTERS
        report_nsaa_sets(neutral_nsaa_flux,nsaa_flux,os);
#endif
      }

      {
        // Report and crudely analyze the stationary distro for each model
        // category. Within each model category, analyze each branch category
        // set
        //
        const rates_func_options opt(bopt,bcs);

        std::ostringstream oss;
        oss << "Cat: " << cm().cat_label(c) << " Branch set: " << branch_set_label;
        const std::string cat_type(oss.str());

        prob_t stat_nscodon[NSCODON::SIZE];
        std::fill(stat_nscodon,stat_nscodon+NSCODON::SIZE,0.);

        stat_pdistro_nscodon(stat_nscodon,opt);

        print_section_divider(os);

        {
          prob_t stat_nuc[NUC::SIZE];
          nscodon_pdf_2_nuc_pdf(stat_nscodon,stat_nuc);

          os << "Nucleotide stationary distribution for " << cat_type << " ::\n\n";

          os << "NUST NUC:    ";
          for(unsigned i(0);i<NUC::SIZE;++i) {
            for(unsigned j(0);j<((pdistro_row_write_width+2)-1);++j) os << " ";
            os << NUC::syms[i];
          }
          os << "\n";

          os << "NUST STAT:   ";
          pdistro_row_write(stat_nuc,NUC::SIZE,os,pdistro_row_write_width+2,5);

          os << "NUST OBS:    ";
          pdistro_row_write(bg_pdistro_nuc_cat(c),NUC::SIZE,os,pdistro_row_write_width+2,5);

          os << std::setprecision(7);
          const smlfloat stdot(normalized_dot(stat_nuc,bg_pdistro_nuc_cat(c),NUC::SIZE));
          os << "NUST dot(STAT,OBS):  " << stdot << "\n";

          const smlfloat stre(rel_ent(stat_nuc,bg_pdistro_nuc_cat(c),NUC::SIZE));
          os << "NUST re(STAT,OBS):   " << stre << "\n";
          os << "\n";
        }

        {
          prob_t stat_nsaa[NSAA::SIZE];
          nscodon_pdf_2_nsaa_pdf(stat_nscodon,stat_nsaa);

          os << "Amino acid stationary distribution for " << cat_type << " ::\n\n";

          os << "AAST AA:    ";
          for(unsigned i(0);i<NSAA::SIZE;++i) {
            for(unsigned j(0);j<(pdistro_row_write_width-1);++j) os << " ";
            os << NSAA::syms[i];
          }
          os << "\n";

          os << "AAST STAT:  ";
          pdistro_row_write(stat_nsaa,NSAA::SIZE,os);

          os << "AAST OBS:   ";
          pdistro_row_write(bg_pdistro_nsaa_cat,NSAA::SIZE,os);

          os << std::setprecision(7);
          const smlfloat stdot(normalized_dot(stat_nsaa,bg_pdistro_nsaa_cat,NSAA::SIZE));
          os << "AAST dot(STAT,OBS):  " << stdot << "\n";

          const smlfloat stre(rel_ent(stat_nsaa,bg_pdistro_nsaa_cat,NSAA::SIZE));
          os << "AAST re(STAT,OBS):   " << stre << "\n";
          os << "\n";
        }

        {
          os << "Codon stationary distribution for " << cat_type << " ::\n\n";

          os << std::setprecision(7);
          const smlfloat stdot(normalized_dot(stat_nscodon,bg_pdistro_nscodon_cat(c),NSCODON::SIZE));

          os << "NCST dot(STAT,OBS):  " << stdot << "\n";

          const smlfloat stre(rel_ent(stat_nscodon,bg_pdistro_nscodon_cat(c),NSCODON::SIZE));
          os << "NCST re(STAT,OBS):   " << stre << "\n";
          os << "\n";
        }
      }
    }
  }
}

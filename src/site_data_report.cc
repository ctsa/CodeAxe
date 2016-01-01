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

// $Id: site_data_report.cc 1195 2008-04-01 02:46:43Z ctsa $

/// \file

#include "simple_util.h"
#include "site_data.h"
#include "site_data_report.h"
#include "site_data_stats.h"
#include "stationary_pdistro.h"
#include "util/bio/bioseq_util.h"
#include "util/bio/bioseq_util_pdf_conversions.h"

#include <ostream>
#include <sstream>

using namespace std;



#if 0
static
void
report_nscodon_distro(const prob_t nsp[],
                      const char* label,
                      std::ostream& os){

  os << label << "\n";
  os << "codon distro:\n";

  NSCODON::print nscodon_syms[NSCODON::SIZE];
  for(unsigned i(0);i<NSCODON::SIZE;++i){ nscodon_syms[i] = NSCODON::print(static_cast<NSCODON::index_t>(i)); }

  pdistro_report(nsp,NSCODON::SIZE,nscodon_syms,os);
}
#endif



static
void
dump_sites(const site_data& sd,
           std::ostream& os){

  const unsigned n_taxa(sd.n_taxa());
  const SITE_MODEL::index_t sm(sd.sm);

  for(unsigned t(0);t<n_taxa;++t) {
    os << sd.taxid.getstr(t) << " ";
  }
  os << "\n";

  site_count_map::const_iterator i,i_end;
  i=sd.site_count.begin();
  i_end=sd.site_count.end();

  for(;i!=i_end;++i){
    site_count_data_type::const_iterator j,j_end;
    j=i->second.begin();
    j_end=i->second.end();

    unsigned count(0);
    for(;j!=j_end;++j){ count += j->second; }

    for(unsigned t(0);t<n_taxa;++t){
      const unsigned s(i->first.get_taxid(t));
      os << SITE_MODEL::print(sm,s,true) << " ";
    }

    os << count << "\n";
  }
}



template <typename SYM_TYPE>
void
report_distro_set(const unsigned n_sets,
                  const unsigned n_states,
                  const prob_t* pdf, // [n_groups*n_states]
                  const name_id_lup& labels,
                  const SYM_TYPE* syms,
                  std::ostream& os){

  for(unsigned i(0); i<n_sets;++i){
    if(i!=0) os << " ";
    os << labels.getstr(i);
  }
  os << "\n";

  os << std::setprecision(5) << std::left << fixed;
  for(unsigned i(0);i<n_states;++i){
    os << syms[i] << " : ";
    for(unsigned j(0);j<n_sets;++j){
      os << std::setw(8) << pdf[i+j*n_states] << " ";
    }
    os << "\n";
  }
  os << "\n";
}



static
void
site_data_synonymous_site_count(const site_data& sd,
                                bool& is_coding,
                                unsigned& synonymous_site_count){

  synonymous_site_count = 0;
  is_coding = SITE_MODEL::is_coding(sd.sm);
  if(! is_coding) return;

  const unsigned n_seqs(sd.n_taxa());

  site_count_map::const_iterator i,i_end=sd.site_count.end();

  for(i=sd.site_count.begin();i!=i_end;++i){
    site_count_data_type::const_iterator j,j_end=i->second.end();

    // check if this is a synon site:
    bool is_synon(true);
    {
      NSAA::index_t last_aa(NSAA::X);
      for(unsigned t(0);t<n_seqs;++t){
        const unsigned state(i->first.get_taxid(t));
        const NSCODON::index_t nsc(decode_nscodon(sd.sm,state));
        const NSAA::index_t aa(codon_trans_known(nsc));
        if(t!=0 && aa!=last_aa) {
          is_synon=false;
          break;
        }
        last_aa=aa;
      }
    }

    // sum over all categories:
    if(is_synon){
      for(j=i->second.begin();j!=j_end;++j) synonymous_site_count += j->second;
    }
  }
}


#include <set>


static
unsigned
site_data_group_count(const site_data& sd){

  std::set<unsigned> active_groups;

  site_count_map::const_iterator i,i_end=sd.site_count.end();

  for(i=sd.site_count.begin();i!=i_end;++i){
    const site_code sc(i->first);
    site_count_data_type::const_iterator j,j_end=i->second.end();

    for(j=i->second.begin();j!=j_end;++j){
      const unsigned group_id(static_cast<unsigned>(j->first.group_id));

      active_groups.insert(group_id);
    }
  }

  return active_groups.size();
}



static
void
site_data_write_counts(const site_data& sd,
                       std::ostream& os){

  unsigned site_count(0);
  unsigned conserved_site_count(0);
  unsigned unique_site_count(0);
  site_data_counts(sd,site_count,conserved_site_count,unique_site_count);

  os << "total sites: " << site_count << "\n";
  os << "conserved sites: " << conserved_site_count << "\n";
  os << "unique sites: " << unique_site_count << "\n";

  bool is_coding(false);
  unsigned synonymous_site_count(0);
  site_data_synonymous_site_count(sd,is_coding,synonymous_site_count);
  if(is_coding) os << "synonymous sites: " << synonymous_site_count << "\n";
  os << "\n";

  const unsigned group_count(site_data_group_count(sd));
  os << "total groups: " << group_count << "\n\n";

  const unsigned n_classes(sd.data_class_size());
  std::vector<unsigned> data_class_count(n_classes,0);
  site_data_2_data_class_count(sd,n_classes,data_class_count);
  os << "total sites per data class:\n";
  for(unsigned c(0);c<n_classes;++c) {
    os << sd.data_class_labels.getstr(c) << " " << data_class_count[c] << "\n";
  }
  os << "\n\n";
}



static
void
report_nscodon_pdf_set(const unsigned n_sets,
                       const prob_t* nscodon_pdf_set,
                       const name_id_lup& labels,
                       std::ostream& os){

  { // nucs:
    prob_t* nuc_pdf_set = new prob_t[n_sets*NUC::SIZE];
    for(unsigned i(0);i<n_sets;++i) nscodon_pdf_2_nuc_pdf(nscodon_pdf_set+i*NSCODON::SIZE,nuc_pdf_set+i*NUC::SIZE);

    report_distro_set<char>(n_sets,NUC::SIZE,nuc_pdf_set,labels,NUC::syms,os);

    // nucs by codon position:
    for(unsigned k(0);k<CODON::BASE_SIZE;++k){
      for(unsigned i(0);i<n_sets;++i) nscodon_pdf_2_nuc_pdf_pos(nscodon_pdf_set+i*NSCODON::SIZE,k,nuc_pdf_set+i*NUC::SIZE);
      os << "CODON POSITION: " << k << "\n";
      report_distro_set<char>(n_sets,NUC::SIZE,nuc_pdf_set,labels,NUC::syms,os);
    }
    delete [] nuc_pdf_set;
  }

  // aas
  prob_t* nsaa_pdf_set = new prob_t[n_sets*NSAA::SIZE];
  for(unsigned i(0);i<n_sets;++i) nscodon_pdf_2_nsaa_pdf(nscodon_pdf_set+i*NSCODON::SIZE,nsaa_pdf_set+i*NSAA::SIZE);

  report_distro_set<char>(n_sets,NSAA::SIZE,nsaa_pdf_set,labels,NSAA::syms,os);
  delete [] nsaa_pdf_set;

  // codons
  NSCODON::print nscodon_syms[NSCODON::SIZE];
  for(unsigned i(0);i<NSCODON::SIZE;++i){ nscodon_syms[i] = NSCODON::print(static_cast<NSCODON::index_t>(i),true); }

  report_distro_set<NSCODON::print>(n_sets,NSCODON::SIZE,nscodon_pdf_set,labels,nscodon_syms,os);
}




static
void
analyze_aa_from_nscodon_pdf(const smlfloat* nscodon_pdf,
                            std::ostream& os){

  static const smlfloat nsexp(1./static_cast<smlfloat>(NSCODON::SIZE));
  static bool is_init(false);
  static map<NSAA::index_t,unsigned> nsaa_codon_count;

  // get codons per aa:
  if(! is_init){
    for(unsigned i(0);i<NSCODON::SIZE;++i){
      NSCODON::index_t ni = static_cast<NSCODON::index_t>(i);
      nsaa_codon_count[codon_trans_known(ni)]++;
    }
    is_init=true;
  }

  // get aa distro
  smlfloat nsaa_pdf[NSAA::SIZE];
  nscodon_pdf_2_nsaa_pdf(nscodon_pdf,nsaa_pdf);

  // get neutral expected aa distro from nuc bg
  smlfloat nuc_pdf[NSAA::SIZE];
  nscodon_pdf_2_nuc_pdf(nscodon_pdf,nuc_pdf);

  smlfloat nscodon_expect_pdf[NSCODON::SIZE];
  nuc_pdf_2_nscodon_pdf(nuc_pdf,nscodon_expect_pdf);

  smlfloat nsaa_expect_pdf[NSCODON::SIZE];
  nscodon_pdf_2_nsaa_pdf(nscodon_expect_pdf,nsaa_expect_pdf);


  os << "AA    OBS    EXP    ln2(O/E)    EXP_NUC_BG    ln2(O/E_N)\n";
  for(unsigned i(0);i<NSAA::SIZE;++i){
    // get naive (uniform nuc distro) neutral aa expectation:
    const smlfloat nsaa_exp = nsexp*nsaa_codon_count[static_cast<NSAA::index_t>(i)];
    os << NSAA::syms[i] << " : " << nsaa_pdf[i] << "   ";

    os << nsaa_exp << " ";
    os.setf(ios::showpos);
    os << std::log(nsaa_pdf[i]/nsaa_exp)/M_LN2 << "   ";
    os.unsetf(ios::showpos);

    os << nsaa_expect_pdf[i] << " ";
    os.setf(ios::showpos);
    os << std::log(nsaa_pdf[i]/nsaa_expect_pdf[i])/M_LN2 << "\n";
    os.unsetf(ios::showpos);
  }
  os << "\n";
}


static
void
state_pdistro_to_nscodon_pdistro(const SITE_MODEL::index_t sm,
                                 const prob_t* state_pdistro,
                                 prob_t* nscodon_pdistro) {

  const BIO_PDISTRO::index_t i(BIO_PDISTRO::convert_from_site_model(sm));
  bio_pdistro_convert(i,BIO_PDISTRO::NSCODON,state_pdistro,nscodon_pdistro);
}


static
void
analyze_nscodon_pdf(const smlfloat* nscodon_pdf,
                    std::ostream& os){

  static const smlfloat nsexp = 1./static_cast<smlfloat>(NSCODON::SIZE);

  smlfloat nuc_pdf[NSAA::SIZE];
  nscodon_pdf_2_nuc_pdf(nscodon_pdf,nuc_pdf);

  smlfloat nscodon_expect_pdf[NSCODON::SIZE];
  nuc_pdf_2_nscodon_pdf(nuc_pdf,nscodon_expect_pdf);

  os << "CODON    OBS      EXP    ln2(O/E)      EXP_NUC_BG    ln2(O/E_N)\n";
  for(unsigned i(0);i<NSCODON::SIZE;++i){
    NSCODON::index_t ni = static_cast<NSCODON::index_t>(i);
    os << NSCODON::print(ni,true) << " : " << nscodon_pdf[i] << "   "
       << nsexp << " ";
    os.setf(ios::showpos);
    os << std::log(nscodon_pdf[i]/nsexp)/M_LN2 << "   ";
    os.unsetf(ios::showpos);
    os << nscodon_expect_pdf[i] << " ";
    os.setf(ios::showpos);
    os << std::log(nscodon_pdf[i]/nscodon_expect_pdf[i])/M_LN2 << "\n";
    os.unsetf(ios::showpos);
  }
  os << "\n";
}



void
site_data_report(const site_data& sd,
                 const bool is_pretty_print_data,
                 std::ostream& os){

  const SITE_MODEL::index_t sm(sd.sm);

  os << "site-type: " << SITE_MODEL::label(sm) << "\n";
  os << "\n";

  site_data_write_counts(sd,os);

  const bool is_c4_site(sm == SITE_MODEL::NSC4PRE ||
                        sm == SITE_MODEL::NSC4POST);

  const unsigned n_org(sd.n_taxa());
  const unsigned n_states(SITE_MODEL::state_size(sm));

  if(sd.sm == SITE_MODEL::NUC){
    prob_t* nuc_pdf = new prob_t[n_org*NUC::SIZE];
    site_data_2_taxa_state_pdf(sd,NUC::SIZE,nuc_pdf);
    report_distro_set<char>(n_org,NUC::SIZE,nuc_pdf,sd.taxid,NUC::syms,os);
    delete [] nuc_pdf;
    return;

  } else if(sd.sm == SITE_MODEL::DINUC){
    DINUC::print dinuc_syms[DINUC::SIZE];
    for(unsigned i(0);i<DINUC::SIZE;++i){ dinuc_syms[i] = DINUC::print(i); }
    {
      prob_t* dinuc_pdf = new prob_t[n_org*DINUC::SIZE];
      site_data_2_taxa_state_pdf(sd,DINUC::SIZE,dinuc_pdf);
      report_distro_set<DINUC::print>(n_org,DINUC::SIZE,dinuc_pdf,sd.taxid,dinuc_syms,os);
      delete [] dinuc_pdf;
    }
    { // more detailed combo and margin reports:
      prob_t all_dinuc_pdf[DINUC::SIZE];
      site_data_2_state_pdf(sd,DINUC::SIZE,all_dinuc_pdf);

      os << "Combined taxa pdf:\n";
      pdistro_report(all_dinuc_pdf,DINUC::SIZE,dinuc_syms,os);

      prob_t all_dinuc_cond_pdf[DINUC::SIZE];
      dinuc_pdf_2_conditioned_dinuc_pdf(all_dinuc_pdf,0,all_dinuc_cond_pdf);

      prob_t n0_stationary_pdistro[NUC::SIZE];
      get_stationary_pdistro_from_pmatrix(n0_stationary_pdistro,all_dinuc_cond_pdf,NUC::SIZE);

      os << "nx+n0 : p(n0|nx):\n";
      pdistro_report(all_dinuc_cond_pdf,DINUC::SIZE,dinuc_syms,os);

      dinuc_pdf_2_conditioned_dinuc_pdf(all_dinuc_pdf,1,all_dinuc_cond_pdf);
      os << "nx+n0 : p(nx|n0):nuc \n";
      pdistro_report(all_dinuc_cond_pdf,DINUC::SIZE,dinuc_syms,os);


      prob_t all_nuc_pdf[NUC::SIZE];
      dinuc_pdf_2_nuc_pdf_pos1(all_dinuc_pdf,all_nuc_pdf);
      os << "n0 pdf:\n";
      pdistro_report(all_nuc_pdf,NUC::SIZE,NUC::syms,os);

      os << "n0 stationary pdf:\n";
      pdistro_report(n0_stationary_pdistro,NUC::SIZE,NUC::syms,os);
    }
    return;
  }

  if(! SITE_MODEL::is_coding(sm)) return;

  { // print out cat specific state distros:

    const unsigned n_classes(sd.data_class_size());
    if(n_classes>1){
      simple_array<prob_t> class_state_pdf(n_classes*n_states);

      site_data_2_data_class_state_pdf(sd,n_classes,n_states,class_state_pdf.ptr());

      // print out cat specific codon distros:
      if(sm == SITE_MODEL::NSCODON){
        report_nscodon_pdf_set(n_classes,class_state_pdf.ptr(),sd.data_class_labels,os);
      } else {
        simple_array<prob_t> class_nscodon_pdf(n_classes*NSCODON::SIZE);

        for(unsigned c(0);c<n_classes;++c){
          state_pdistro_to_nscodon_pdistro(sm,
                                           class_state_pdf.ptr()+c*n_states,
                                           class_nscodon_pdf.ptr()+c*NSCODON::SIZE);
        }
        report_nscodon_pdf_set(n_classes,class_nscodon_pdf.ptr(),sd.data_class_labels,os);

        simple_array<SITE_MODEL::print> syms(n_states);
        for(unsigned i(0);i<n_states;++i) syms[i] = SITE_MODEL::print(sd.sm,i);

        report_distro_set<SITE_MODEL::print>(n_classes,n_states,class_state_pdf.ptr(),sd.data_class_labels,syms.ptr(),os);
      }
    }
  }


  { // print out org specific state distros:
    simple_array<prob_t> org_state_pdf(n_org*n_states);

    site_data_2_taxa_state_pdf(sd,n_states,org_state_pdf.ptr());

    // print out org specific codon distros:
    if(sm == SITE_MODEL::NSCODON){
      report_nscodon_pdf_set(n_org,org_state_pdf.ptr(),sd.taxid,os);
    } else {
      simple_array<prob_t> org_nscodon_pdf(n_org*NSCODON::SIZE);

      for(unsigned t(0);t<n_org;++t){
        state_pdistro_to_nscodon_pdistro(sm,
                                         org_state_pdf.ptr()+t*n_states,
                                         org_nscodon_pdf.ptr()+t*NSCODON::SIZE);
      }
      report_nscodon_pdf_set(n_org,org_nscodon_pdf.ptr(),sd.taxid,os);


      simple_array<SITE_MODEL::print> syms(n_states);
      for(unsigned i(0);i<n_states;++i) syms[i] = SITE_MODEL::print(sd.sm,i);

      report_distro_set<SITE_MODEL::print>(n_org,n_states,org_state_pdf.ptr(),sd.taxid,syms.ptr(),os);
    }
  }

  { // report combined-taxa summary:
    simple_array<prob_t> state_pdf(n_states);

    site_data_2_state_pdf(sd,n_states,state_pdf.ptr());

    prob_t nscodon_pdf[NSCODON::SIZE];
    state_pdistro_to_nscodon_pdistro(sm,
                                     state_pdf.ptr(),
                                     nscodon_pdf);

    analyze_aa_from_nscodon_pdf(nscodon_pdf,os);
    analyze_nscodon_pdf(nscodon_pdf,os);

    if(is_c4_site){
      // supplemental ultrabug analysis
      prob_t nx_pdf[NUC::SIZE];
      prob_t n1_pdf[NUC::SIZE];
      prob_t n2_pdf[NUC::SIZE];
      prob_t nxn2_pdf[NUC::SIZE*NUC::SIZE];

      prob_t n1n2_pdf[NUC::SIZE*NUC::SIZE];
      for(unsigned i(0);i<NUC::SIZE;++i){
        nx_pdf[i] = 0.;
        n1_pdf[i] = 0.;
        n2_pdf[i] = 0.;
      }
      for(unsigned i(0);i<(NUC::SIZE*NUC::SIZE);++i){
        nxn2_pdf[i] = 0.;
        n1n2_pdf[i] = 0.;
      }

      for(unsigned i(0);i<NSC4::SIZE;++i){
        NUC::index_t nx;
        NSCODON::index_t c;
        NSC4::decode(nx,c,i);

        NUC::index_t n[3];
        NSCODON::decode(n,c);
        const prob_t val(state_pdf[i]);
        nx_pdf[nx] += val;
        n1_pdf[n[1]] += val;
        n2_pdf[n[2]] += val;
        nxn2_pdf[n[2]+NUC::SIZE*nx] += val;
        n1n2_pdf[n[2]+NUC::SIZE*n[1]] += val;
      }

      os << "c4 nx distro:\n";
      pdistro_report(nx_pdf,NUC::SIZE,NUC::syms,os);

      os << "c4 distros:\n"
         << "nx-c: obs exp:\n";
      for(unsigned i(0);i<NUC::SIZE;++i){
        for(unsigned j(0);j<NSCODON::SIZE;++j){
          const NSC4::index_t ni(NSC4::encode(static_cast<NUC::index_t>(i),
                                              static_cast<NSCODON::index_t>(j)));
          const smlfloat obs(state_pdf[ni]);
          const smlfloat exp(nx_pdf[i]*nscodon_pdf[j]);
          os << SITE_MODEL::print(sd.sm,ni) << ": "
             << obs << " "
             << exp << " "
             << std::log(obs/exp) << "\n";
        }
      }
      os << "\n";


      os << "c4 n1 distro:\n";
      pdistro_report(n1_pdf,NUC::SIZE,NUC::syms,os);

      os << "c4 n2 distro:\n";
      pdistro_report(n2_pdf,NUC::SIZE,NUC::syms,os);

      {
        DINUC::print dinuc_syms[DINUC::SIZE];
        for(unsigned i(0);i<DINUC::SIZE;++i){ dinuc_syms[i] = DINUC::print(i); }

        prob_t all_dinuc_cond_pdf[DINUC::SIZE];
        dinuc_pdf_2_conditioned_dinuc_pdf(n1n2_pdf,0,all_dinuc_cond_pdf);

        os << "c4: n1+n2 : p(n2|n1):\n";
        pdistro_report(all_dinuc_cond_pdf,DINUC::SIZE,dinuc_syms,os);

        dinuc_pdf_2_conditioned_dinuc_pdf(n1n2_pdf,1,all_dinuc_cond_pdf);

        os << "c4: n1+n2 : p(n1|n2):\n";
        pdistro_report(all_dinuc_cond_pdf,DINUC::SIZE,dinuc_syms,os);
      }
#if 0
      os << "\n"
         << "double nuc distros:\n"
         << "N-M: obs exp:\n";
      for(unsigned i(0);i<NUC::SIZE;++i){
        for(unsigned j(0);j<NUC::SIZE;++j){
          smlfloat obs=nxn2_pdf[j+NUC::SIZE*i];
          smlfloat exp=nx_pdf[i]*n2_pdf[j];
          os << NUC::syms[i] << "-" << NUC::syms[j] << ": "
             << obs << " "
             << exp << " "
             << std::log(obs/exp) << "\n";
        }
      }
#endif
    }
  }

  // dump all data in readable format:
  if(is_pretty_print_data) {
    os << "PRETTY_PRINT_SITES:\n";
    dump_sites(sd,os);
  }
}




#include "cat_info.h"

#include <fstream>


static
void
site_data_2_single_class_site_data(const site_data& sd,
                                   const unsigned n_cats,
                                   const unsigned target_data_class_id,
                                   site_data& sd_out){

  sd_out.clear();

  sd_out.taxid=sd.taxid;
  sd_out.sm=sd.sm;
  sd_out.group_label=sd.group_label;

  const std::string& s(sd.data_class_labels.getstr(target_data_class_id));
  sd_out.data_class_labels.assignid(s);

  site_count_map::const_iterator i,i_end=sd.site_count.end();

  for(i=sd.site_count.begin();i!=i_end;++i){
    const site_code sc(i->first);
    site_count_data_type::const_iterator j,j_end=i->second.end();

    for(j=i->second.begin();j!=j_end;++j){
      const unsigned group_id(static_cast<unsigned>(j->first.group_id));
      const unsigned data_class_id(static_cast<unsigned>(j->first.data_class_id));

      if(data_class_id>=n_cats) continue;

      if(data_class_id!=target_data_class_id) continue;

      const unsigned count(j->second);

      sd_out.site_count[sc][site_data_code(group_id,UNASSIGNED_CAT_ID)] += count;
    }
  }
}


void
site_data_cat_split(const site_data& sd,
                    const char* outtag){


  site_data sd_out;

  const unsigned n_classes(sd.data_class_size());
  for(unsigned i(0);i<n_classes;++i){
    site_data_2_single_class_site_data(sd,n_classes,i,sd_out);

    std::stringstream oss;
    oss << outtag << "_data_class_" << i;

    std::ofstream fos(oss.str().c_str());

    sd_out.store_state(fos);
  }
}

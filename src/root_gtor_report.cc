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

// $Id: root_gtor_report.cc 1210 2008-05-20 23:17:36Z ctsa $

/// \file

#include "cat_manager.h"
#include "root_gtor.h"
#include "stationary_pdistro.h"
#include "util/bio/bioseq_util_io.h"
#include "util/bio/bioseq_util_pdf_conversions.h"

#include <algorithm>
#include <ostream>



/// keep things consistent
static
void
format_root_header(const unsigned root_distro_cat_no,
                   const char* root_label,
                   std::ostream& os){

  os << root_label << " prior prob: root_distro_cat_no " << static_cast<int>(root_distro_cat_no);
}



static
void
report_nuc_distro_from_nscodon(const prob_t nsp[],
                               const char* label,
                               std::ostream& os){

  os << label << "\n";
  os << "nuc distro:\n";

  prob_t nuc_pdf[NUC::SIZE];
  nscodon_pdf_2_nuc_pdf(nsp,nuc_pdf);
  pdistro_report(nuc_pdf,NUC::SIZE,NUC::syms,os);
}



static
void
report_nuc_distro_from_nscodon_pos(const prob_t nscodon_pdf[],
                                   const char* label,
                                   const unsigned pos,
                                   std::ostream& os){

  os << label << "\n";
  os << "nuc distro from codon pos " << pos << ":\n";

  prob_t nuc_pdf[NUC::SIZE];
  nscodon_pdf_2_nuc_pdf_pos(nscodon_pdf,pos,nuc_pdf);
  pdistro_report(nuc_pdf,NUC::SIZE,NUC::syms,os);
}



static
void
report_nsaa_distro_from_nscodon(const prob_t nsp[],
                                const char* label,
                                std::ostream& os){
  os << label << "\n";
  os << "aa distro:\n";

  prob_t nsaa_pdf[NSAA::SIZE];
  nscodon_pdf_2_nsaa_pdf(nsp,nsaa_pdf);
  pdistro_report(nsaa_pdf,NSAA::SIZE,NSAA::syms,os);
}



static
void
report_root_dot_re(const char* tag,
                   const prob_t* root_pdistro,
                   const prob_t* bg_pdistro,
                   const unsigned size,
                   std::ostream& os){

  os << std::setprecision(7);
  os << tag << "RO dot(ROOT,OBS):  ";
  os << normalized_dot(root_pdistro,bg_pdistro,size) << "\n";

  os << tag << "RO re(ROOT,OBS):   ";
  os << rel_ent(root_pdistro,bg_pdistro,size) << "\n";
  os << "\n";
}



const unsigned pdistro_row_write_width(7);
const unsigned pdistro_row_write_prec(3);


static
void
pdistro_row_write(const prob_t* pdistro,
                  const unsigned N,
                  std::ostream& os,
                  const unsigned width = pdistro_row_write_width,
                  const unsigned prec = pdistro_row_write_prec){

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



static
void
report_compare_sec(const char* distro_label,
                   const char* site_label,
                   const unsigned wadjust,
                   const char* syms,
                   const char* tag,
                   const prob_t* root_pdistro,
                   const prob_t* bg_pdistro,
                   const unsigned size,
                   std::ostream& os){

  const unsigned new_width(pdistro_row_write_width+wadjust);
  const unsigned new_prec(pdistro_row_write_prec+wadjust);

  os << distro_label << "\n";
  os << "root comparison to obs distro (" << site_label << "):\n";

  os << tag << "RO:      ";
  for(unsigned k(0);k<size;++k) {
    for(unsigned j(0);j<(new_width-1);++j) os << " ";
    os << syms[k];
  }
  os << "\n";

  os << tag << "RO ROOT: ";
  pdistro_row_write(root_pdistro,size,os,new_width,new_prec);

  os << tag << "RO OBS:  ";
  pdistro_row_write(bg_pdistro,size,os,new_width,new_prec);

  report_root_dot_re(tag,root_pdistro,bg_pdistro,size,os);
}



/// hack for c4+ mode only
static
void
print_seq_stat_root(const prob_t* root_pdf,
                    const unsigned ss,
                    std::ostream& os){

  // get sequence stationary root
  //
  static const unsigned repeat_pos(3);
  static const unsigned dep_pos(0);

  static const SITE_MODEL::index_t sm(SITE_MODEL::NSC4POST);

  // get stationary dep nuc distro
  prob_t dep_nuc_stat[NUC::SIZE];
  {
    // extract dependent nuc distro
    prob_t dep_nuc_p[NUC::SIZE*NUC::SIZE];

    std::fill(dep_nuc_p,dep_nuc_p+NUC::SIZE*NUC::SIZE,0.);

    NUC::index_t sn[SITE_MODEL::MAX_BASE_SIZE];
    for(unsigned j(0);j<ss;++j){
      SITE_MODEL::decode_nuc(sm,j,sn);
      dep_nuc_p[sn[repeat_pos]*NUC::SIZE+sn[dep_pos]] += root_pdf[j];
    }

    for(unsigned j(0);j<NUC::SIZE;++j){
      pdistro_norm(dep_nuc_p+j*NUC::SIZE,dep_nuc_p+(j+1)*NUC::SIZE);
    }
    get_stationary_pdistro_from_pmatrix(dep_nuc_stat,dep_nuc_p,NUC::SIZE);
  }

  // get sequence-stationary codon distro
  prob_t sstat_nscodon_pdf[NSCODON::SIZE];
  {
    prob_t nscodon_by_nx_base[NSC4::SIZE];
    for(unsigned j(0);j<NSC4::SIZE;++j) nscodon_by_nx_base[j] = root_pdf[j];

    prob_t* nscodon_by_nx[NUC::SIZE];
    for(unsigned j(0);j<NUC::SIZE;++j){
      nscodon_by_nx[j] = nscodon_by_nx_base+NSC4::get_nscodon_offset(static_cast<NUC::index_t>(j));
    }

    for(unsigned j(0);j<NUC::SIZE;++j){
      pdistro_norm(nscodon_by_nx[j],nscodon_by_nx[j]+NSCODON::SIZE);
    }

    std::fill(sstat_nscodon_pdf,sstat_nscodon_pdf+NSCODON::SIZE,0.);

    for(unsigned k(0);k<NUC::SIZE;++k){
      for(unsigned j(0);j<NSCODON::SIZE;++j){
        sstat_nscodon_pdf[j] += nscodon_by_nx[k][j]*dep_nuc_stat[k];
      }
    }
  }

  print_site_model_distro(sstat_nscodon_pdf,SITE_MODEL::NSCODON,os);
}



void
root_gtor::
report(const char* label,
       std::ostream& os) const {

  const RATE_GTOR_MODEL::index_t ragm(_smls.rate_gtor_model());
  const SITE_MODEL::index_t sm(RATE_GTOR_MODEL::convert_to_site_model(ragm));
  const BIO_PDISTRO::index_t bp(BIO_PDISTRO::convert_from_site_model(sm));

  /// \todo move all this stuff to heap?
  prob_t obs_nscodon_pdistro[NSCODON::SIZE];
  prob_t obs_nsaa_pdistro[NSAA::SIZE];
  prob_t obs_nuc_pdistro[NUC::SIZE];
  prob_t stack_root_nscodon_pdistro[NSCODON::SIZE];
  prob_t stack_root_nuc_pdistro[NUC::SIZE];

  const unsigned ss(state_size());
  const unsigned n_cats(cm().cat_size());

  os << "root_gtor_report:\n\n";

  const unsigned ds(distro_cat_size());
  for(unsigned i(0);i<ds;++i){
    std::string distro_label;
    {
      std::ostringstream oss;
      format_root_header(i,label,oss);
      distro_label = oss.str();
    }

    const prob_t* root_pdistro(distro_cat_state_pdistro(i));

    /// average all obs distros associated with this root distro cat:
    std::vector<prob_t> obs_pdistro(ss,0.);
    {
      std::vector<prob_t> cat_pdistro(n_cats);
      cm().cat_pdistro(cat_pdistro.begin());
      for(unsigned j(0);j<n_cats;++j){
        if(get_distro_cat_no_from_cat_no(j)!=i) continue;
        const prob_t* ocp(_smls.obs_state_distro_cat(j));
        const prob_t cp(cat_pdistro[j]);
        for(unsigned k(0);k<ss;++k) obs_pdistro[k] += ocp[k]*cp;
      }
      pdistro_norm(obs_pdistro.begin(),obs_pdistro.end());
    }

    if(! SITE_MODEL::is_coding(sm)){ // nuc version
      const prob_t* root_nuc_pdistro(root_pdistro);
      if(bp != BIO_PDISTRO::NUC){
        os << distro_label << "\n";
        print_site_model_distro(root_pdistro,sm,os);

        bio_pdistro_convert(bp,BIO_PDISTRO::NUC,root_pdistro,stack_root_nuc_pdistro);
        root_nuc_pdistro=stack_root_nuc_pdistro;
      }

      os << distro_label << "\n";
      print_site_model_distro(root_nuc_pdistro,SITE_MODEL::NUC,os);

      // report differences between root distribution and obs:
      bio_pdistro_convert(bp,BIO_PDISTRO::NUC,obs_pdistro.begin(),obs_nuc_pdistro);
      report_compare_sec(distro_label.c_str(),"nuc",2,NUC::syms,
                         "NU",root_nuc_pdistro,obs_nuc_pdistro,NUC::SIZE,os);

      continue;
    }

    // nscodon version: need to resolve this w/ nuc version!
    bio_pdistro_convert(bp,BIO_PDISTRO::NSCODON,obs_pdistro.begin(),obs_nscodon_pdistro);
    nscodon_pdf_2_nsaa_pdf(obs_nscodon_pdistro,obs_nsaa_pdistro);
    nscodon_pdf_2_nuc_pdf(obs_nscodon_pdistro,obs_nuc_pdistro);

    const prob_t* root_nscodon_pdistro(root_pdistro);

    if(bp != BIO_PDISTRO::NSCODON){
      bio_pdistro_convert(bp,BIO_PDISTRO::NSCODON,root_pdistro,stack_root_nscodon_pdistro);
      root_nscodon_pdistro=stack_root_nscodon_pdistro;
    }

    // nuc
    report_nuc_distro_from_nscodon(root_nscodon_pdistro,distro_label.c_str(),os);

    // nuc by codon pos
    for(unsigned k(0);k<CODON::BASE_SIZE;++k){
      report_nuc_distro_from_nscodon_pos(root_nscodon_pdistro,distro_label.c_str(),k,os);
    }

    // aa
    report_nsaa_distro_from_nscodon(root_nscodon_pdistro,distro_label.c_str(),os);

    // codons
    os << distro_label << "\n";
    print_site_model_distro(root_nscodon_pdistro,SITE_MODEL::NSCODON,os);

    if(sm==SITE_MODEL::NSC4POST) {
      os << distro_label << " (sequence-stationary)\n";
      print_seq_stat_root(root_pdistro,ss,os);
    }

    // report differences between root distribution and obs:
    //
    {
      prob_t root_nuc_pdistro[NUC::SIZE];
      nscodon_pdf_2_nuc_pdf(root_nscodon_pdistro,root_nuc_pdistro);
      report_compare_sec(distro_label.c_str(),"nuc",2,NUC::syms,
                         "NU",root_nuc_pdistro,obs_nuc_pdistro,NUC::SIZE,os);

      prob_t root_nsaa_pdistro[NSAA::SIZE];
      nscodon_pdf_2_nsaa_pdf(root_nscodon_pdistro,root_nsaa_pdistro);
      report_compare_sec(distro_label.c_str(),"nsaa",0,NSAA::syms,
                         "AA",root_nsaa_pdistro,obs_nsaa_pdistro,NSAA::SIZE,os);

      os << distro_label << "\n";
      os << "root comparison to obs distro (nscodon):\n";

      report_root_dot_re("NC",root_nscodon_pdistro,obs_nscodon_pdistro,NSCODON::SIZE,os);
    }
  }
}

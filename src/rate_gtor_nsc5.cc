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

// $Id: rate_gtor_nsc5.cc 1107 2008-01-25 02:40:20Z ctsa $

/// \file

#include "rate_gtor_nsc5.h"
#include "rate_gtor_nscodon_base_util.h"
#include "subs_ml_ptol.h"
#include "util/bio/bioseq_util_pdf_conversions.h"
#include "util/general/io_util.h"

#include <iostream>
#include <string>

#ifdef DEBUG
#define NSC5_DEBUG
#endif

void
rate_gtor_nsc5::
store_state(std::ostream& os) const {

  base_t::store_state(os);

  os << "rate_gtor_nsc5 c5_approx_model " << _c5a << "\n";
}


void
rate_gtor_nsc5::
load_state_internal(std::istream& is) {

  base_t::load_state_internal(is);

  std::string dummy;
  is >> dummy >> dummy;
  _c5a = static_cast<C5_APPROX::index_t>(read_int(is));
}


void
rate_gtor_nsc5::
rates(smlfloat rts[],
      const rates_func_options& opt) const {
  if( opt.is_use_submodel && (_c5a==C5_APPROX::C4_JOINT || _c5a==C5_APPROX::C4_GMEAN)){

    SITE_MODEL::index_t sm;
    if(opt.submodel_no == 0){ sm=SITE_MODEL::NSC4PRE; }
    else                    { sm=SITE_MODEL::NSC4POST; }

    rates_nscodon_context(*this,sm,rts,opt);
  } else {
    rates_nscodon_context(*this,rts,opt);
  }
}


void
rate_gtor_nsc5::
rates_variance(irv_t<smlfloat> rts[],
               const rates_func_options& opt) const {
  if( opt.is_use_submodel && (_c5a==C5_APPROX::C4_JOINT || _c5a==C5_APPROX::C4_GMEAN)){

    SITE_MODEL::index_t sm;
    if(opt.submodel_no == 0){ sm=SITE_MODEL::NSC4PRE; }
    else                    { sm=SITE_MODEL::NSC4POST; }

    rates_nscodon_context(*this,sm,rts,opt);
  } else {
    rates_nscodon_context(*this,rts,opt);
  }
}


void
rate_gtor_nsc5::
submodel_pdistro_reduction(const prob_t* in,
                           unsigned s,
                           prob_t* out) const {
  if(_c5a==C5_APPROX::C4_JOINT || _c5a==C5_APPROX::C4_GMEAN){
    const bool is_c4_pre(s==0);
    nsc5_pdf_2_nsc4_pdf(in,is_c4_pre,out);
  } else {
    base_t::submodel_pdistro_reduction(in,s,out);
  }
}


void
rate_gtor_nsc5::
submodel_fuse_pdistros(const prob_t* const * in,
                       prob_t* out) const {
  if       (_c5a==C5_APPROX::C4_GMEAN){

    const prob_t* pre_c4(in[0]);
    const prob_t* post_c4(in[1]);

#ifdef NSC5_DEBUG
    pdistro_check(pre_c4,NSC4::SIZE,SUBS_ML_PTOL);
    pdistro_check(post_c4,NSC4::SIZE,SUBS_ML_PTOL);
#endif

    // get averaged codon distro:
    prob_t acodon[NSCODON::SIZE];
    {
      prob_t post_codon[NSCODON::SIZE];
      nsc4_pdf_2_nscodon_pdf(pre_c4,acodon);
      nsc4_pdf_2_nscodon_pdf(post_c4,post_codon);

      for(unsigned i(0);i<NSCODON::SIZE;++i){
        acodon[i] = (acodon[i]+post_codon[i])/2.;
      }
    }

    //
    for(unsigned c(0);c<NSCODON::SIZE;++c){
      const NSCODON::index_t ci(static_cast<NSCODON::index_t>(c));

      prob_t n0sum(0),n4sum(0);
      for(unsigned n(0);n<NUC::SIZE;++n){
        NSC4::index_t c4i(NSC4::encode(static_cast<NUC::index_t>(n),ci));
        n0sum += pre_c4[c4i];
        n4sum += post_c4[c4i];
      }
      if(std::abs(n0sum)<=0) n0sum=1.;
      if(std::abs(n4sum)<=0) n4sum=1.;

      for(unsigned n0(0);n0<NUC::SIZE;++n0){
        const NUC::index_t n0i(static_cast<NUC::index_t>(n0));
        NSC4::index_t pre4i(NSC4::encode(n0i,ci));
        for(unsigned n4(0);n4<NUC::SIZE;++n4){
        const NUC::index_t n4i(static_cast<NUC::index_t>(n4));
          NSC4::index_t post4i(NSC4::encode(n4i,ci));
          NSC5::index_t c5i(NSC5::encode(n0i,n4i,ci));
          out[c5i] = acodon[ci]*(pre_c4[pre4i]/n0sum)*(post_c4[post4i]/n4sum);
        }
      }
    }

#ifdef NSC5_DEBUG
    pdistro_check(out,NSC5::SIZE,SUBS_ML_PTOL);
#endif

  } else if(_c5a==C5_APPROX::C4_JOINT){
    die("no fuse options for c4-point approximation");
  } else {
    base_t::submodel_fuse_pdistros(in,out);
  }
}


void
rate_gtor_nsc5::
submodel_state_reduction_map(unsigned s,
                             unsigned* out) const {
  if(_c5a==C5_APPROX::C4_JOINT || _c5a==C5_APPROX::C4_GMEAN){
    const bool is_c4_pre(s==0);
    NUC::index_t n_pre,n_post,nx;
    NSCODON::index_t c;
    for(unsigned i(0);i<NSC5::SIZE;++i){
      NSC5::decode(n_pre,n_post,c,static_cast<NSC5::index_t>(i));
      if(is_c4_pre) nx = n_pre;
      else          nx = n_post;
      out[i] = NSC4::encode(nx,c);
    }
    // ambiguous state remaps to ambiguous sub-state
    out[NSC5::SIZE] = NSC4::SIZE;
  } else {
    base_t::submodel_state_reduction_map(s,out);
  }
}

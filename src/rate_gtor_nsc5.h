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

// $Id: rate_gtor_nsc5.h 1133 2008-01-29 22:42:35Z ctsa $

/// \file

#ifndef __RATE_GTOR_NSC5_H
#define __RATE_GTOR_NSC5_H


#include "condition_func_nsc5.h"
#include "condition_func_nsc5_center.h"
#include "condition_func_nsc4_joint1.h"
#include "condition_func_nsc4_joint2.h"
#include "condition_func_nsc4.h"
#include "rate_gtor_nscodon_base.h"

#include "util/math/array_util.h"


struct rate_gtor_nsc5 : public rate_gtor_nscodon_base {

  typedef rate_gtor_nscodon_base base_t;
  typedef rate_gtor_nsc5 self_t;

  rate_gtor_nsc5(const rate_gtor_options& ropt,
                 const rate_gtor_nuc_options& nopt,
                 const rate_gtor_nscodon_options& copt,
                 const C5_APPROX::index_t c5a,
                 const rate_gtor_sml_share& smls)
    : base_t(ropt,nopt,copt,smls), _c5a(c5a) {

    if(_c5a !=C5_APPROX::NONE && model_type() != RATE_GTOR_MODEL::C5)
      die("c5 approx mode only valid for centered c5 site model");
  }

  rate_gtor_nsc5(const self_t& s,
                 const rate_gtor_sml_share& smls)
    : base_t(s,smls), _c5a(s._c5a) {}

  ////////////////////////////////////////////////////////////////////////////
private:

  virtual
  rate_gtor* clone(const rate_gtor_sml_share& smls) const { return new self_t(*this,smls); }

  virtual
  unsigned state_size_conditioned() const { return NSCODON::SIZE; }

  virtual
  void rates(smlfloat rates[],
             const rates_func_options& opt) const;

  virtual
  void rates_variance(irv_t<smlfloat> rates[],
                      const rates_func_options& opt) const;

  virtual
  condition_func* condition_func_factory() const {
    if(model_type()==RATE_GTOR_MODEL::C5){
      return new condition_func_nsc5_center();
    } else {
      return new condition_func_nsc5();
    }
  }

  virtual
  void store_state(std::ostream& os) const;

  virtual
  void load_state_internal(std::istream& is);

  // submodel functions enable joint c4 approximation models:
  virtual
  unsigned submodel_size() const {
    if(_c5a==C5_APPROX::C4_JOINT || _c5a==C5_APPROX::C4_GMEAN) return 2;
    else return base_t::submodel_size();
  }

  virtual
  unsigned submodel_state_size(unsigned s) const {
    if(s>=submodel_size()) abort();
    if(_c5a==C5_APPROX::C4_JOINT || _c5a==C5_APPROX::C4_GMEAN) return NSC4::SIZE;
    else return base_t::submodel_state_size(s);
  }

  virtual
  void submodel_adjust_p(prob_t* p,const unsigned n) const {
    if(_c5a==C5_APPROX::C4_GMEAN) array_sqrt(p,n);
    else                          base_t::submodel_adjust_p(p,n);
  }


  virtual
  condition_func* submodel_condition_func_factory(unsigned s) const {
    if       (_c5a==C5_APPROX::C4_JOINT){
      const bool is_c4_pre(s==0);
      if(is_c4_pre){ return new condition_func_nsc4_joint1; }
      else         { return new condition_func_nsc4_joint2; }
    } else if(_c5a==C5_APPROX::C4_GMEAN){
      return new condition_func_nsc4(true);
    } else {
      return base_t::submodel_condition_func_factory(s);
    }
  }

  virtual
  void submodel_pdistro_reduction(const prob_t* in,unsigned s,prob_t* out) const;

  virtual
  void submodel_fuse_pdistros(const prob_t* const * in, prob_t* out) const;

  virtual
  void submodel_state_reduction_map(unsigned s,unsigned* out) const;

  void rates_nsc5(smlfloat rates[],
                  const rates_func_options& opt) const;

  C5_APPROX::index_t _c5a;
};


#endif

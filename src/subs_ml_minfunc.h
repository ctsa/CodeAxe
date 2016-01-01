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

// $Id: subs_ml_minfunc.h 1176 2008-03-26 00:17:56Z ctsa $

/// \file

#ifndef __SUBS_ML_MINFUNC_H
#define __SUBS_ML_MINFUNC_H

#include "simple_util.h"
#include "subs_ml_model.h"
#include "tree_probs.h"
#include "util/general/uncopyable.h"
#include "util/math/minfunc_interface.h"

struct condition_func;
struct estep_dat_t;
struct site_data_fastlup;
struct lhood_model_prep;

extern unsigned min_call_count;
extern unsigned dmin_call_count;



struct subs_ml_minfunc_base : public minfunc_gradient_interface<smlfloat>, private uncopyable {

  subs_ml_minfunc_base(subs_ml_model& m,
                       PARAM_VIEW::index_t pv)
    : _mdl(m), _pv(pv) {}

  virtual
  unsigned dim() const { return po().param_size(_pv); }

  void min_p_to_param(const smlfloat* p) {
    min_p_to_param_x(p,po());
  }

  void param_to_min_p(smlfloat* p) const {
    po().param_state(p,_pv);
  }

  virtual
  smlfloat dval(const smlfloat* p,smlfloat* df);

  virtual
  smlfloat val(const smlfloat* p);

  virtual
  bool is_val_computable(const smlfloat* p);

private:
  virtual
  param_object_base&
  po_from_mdl(subs_ml_model& m) const { return m; }

  virtual
  const param_object_base&
  po_from_mdl(const subs_ml_model& m) const { return m; }

  /// playing games with const:
  ///
  void min_p_to_param_x(const smlfloat* p,
                        param_object_base& po_copy) const {
    check_p(p);
    po_copy.set_param_state(p,_pv);
  }

  virtual bool is_inner_step() const { return false; }

  virtual void get_lnprob(const subs_ml_model& mdl_copy,
                          smlfloat& lnp,
                          smlfloat& lnp_norm) const = 0;

  smlfloat get_deriv_estimate(const smlfloat* p,
                              const smlfloat lnp1_norm,
                              const unsigned param_no,
                              const bool is_plus = true) const;

  void check_p(const smlfloat* p) const;

  bool is_current_val_computable() const;

  param_object_base& po() { return po_from_mdl(mdl()); }
  const param_object_base& po() const { return po_from_mdl(mdl()); }

  subs_ml_model& mdl() { return _mdl; }
  const subs_ml_model& mdl() const { return _mdl; }

  /////////////////////////////////////
  subs_ml_model& _mdl;
  PARAM_VIEW::index_t _pv;
};



struct subs_ml_minfunc : public subs_ml_minfunc_base {

  subs_ml_minfunc(const site_data_fastlup& sdf,
                  subs_ml_model& m,
                  const PARAM_VIEW::index_t pv);

  // defined in subs_ml_minfunc.cc so that auto_ptr to incomplete type works dtors
  // correctly
  ~subs_ml_minfunc();

private:
  virtual
  void get_lnprob(const subs_ml_model& mdl_copy,
                  smlfloat& lnp,
                  smlfloat& lnp_norm) const;

  /////////////////////////////////////////
  const site_data_fastlup& _sdf;
  std::auto_ptr<lhood_model_prep> _nld;
};


#if 0
struct subs_ml_minfunc_em : public subs_ml_minfunc_base {

  subs_ml_minfunc_em(const estep_dat_t& edat,
                     subs_ml_model& m,
                     const PARAM_VIEW::index_t pv)
    : subs_ml_minfunc_base(m,pv), _edat(edat) {}

private:
  virtual
  bool is_inner_step() const { return true; }

  virtual
  void get_lnprob(const subs_ml_model& mdl_copy,
                  smlfloat& lnp,
                  smlfloat& lnp_norm) const;

  /////////////////////////////////////////
  const estep_dat_t& _edat;
};
#endif


/// \brief minfunc for all root categories in root-cycle-mode
///
struct subs_ml_root_minfunc : public subs_ml_minfunc_base {

  subs_ml_root_minfunc(const site_data_fastlup& sdf,
                       subs_ml_model& m,
                       const PARAM_VIEW::index_t pv);

  ~subs_ml_root_minfunc();

  void reset(const subs_ml_model& m);

private:
  virtual
  param_object_base&
  po_from_mdl(subs_ml_model& m) const;

  virtual
  const param_object_base&
  po_from_mdl(const subs_ml_model& m) const;

  virtual
  void get_lnprob(const subs_ml_model& mdl_copy,
                  smlfloat& lnp,
                  smlfloat& lnp_norm) const;

  /////////////////////////////////////////
  bool _is_reset;
  const site_data_fastlup& _sdf;
  std::auto_ptr<lhood_model_prep> _nld;
  simple_matrix3d<smlfloat> _root_spp;
  simple_array<std::auto_ptr<condition_func> > _cf;
};



/// \brief minfunc for a single root category in rootcat-cycle-mode
///
struct subs_ml_rootcat_minfunc : public subs_ml_minfunc_base {

  subs_ml_rootcat_minfunc(const site_data_fastlup& sdf,
                           subs_ml_model& m,
                           const PARAM_VIEW::index_t pv);

  ~subs_ml_rootcat_minfunc();

  void reset(const subs_ml_model& m,
             const unsigned cat_no);

private:
  virtual
  param_object_base&
  po_from_mdl(subs_ml_model& m) const;

  virtual
  const param_object_base&
  po_from_mdl(const subs_ml_model& m) const;

  virtual
  void get_lnprob(const subs_ml_model& mdl_copy,
                  smlfloat& lnp,
                  smlfloat& lnp_norm) const;

  /////////////////////////////////////////
  bool _is_reset;
  unsigned _cat_no;
  const site_data_fastlup& _sdf;
  std::auto_ptr<lhood_model_prep> _nld;
  simple_matrix<smlfloat> _root_spp;
  simple_matrix<smlfloat> _othercat_site_prob;
  std::auto_ptr<condition_func> _cf;
};



#endif

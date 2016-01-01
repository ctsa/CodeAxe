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

// $Id: prob_gtor_obs_derived.h 1118 2008-01-28 20:32:23Z ctsa $

/// \file

#ifndef __PROB_GTOR_OBS_DERIVED_H
#define __PROB_GTOR_OBS_DERIVED_H

#include "prob_gtor_obs.h"


/// \brief simple average between the leaf obs distros
///
struct prob_gtor_obs_state_avg : public prob_gtor_obs {

  typedef prob_gtor_obs_state_avg self_t;
  typedef prob_gtor_obs base_t;

  prob_gtor_obs_state_avg(const prob_gtor_obs_sml_share& smlm,
                          const unsigned cat_no)
    : base_t(smlm,cat_no) {}

  prob_gtor_obs_state_avg(const self_t& s,
                          const prob_gtor_obs_sml_share& smlm)
    : base_t(s,smlm) {}

  virtual self_t* clone(const prob_gtor_obs_sml_share& smlm) const { return new self_t(*this,smlm); }

private:
  virtual void eval_pdistro_update() const {
    _smlm.get_obs_state_avg_root(pdistro_mutable(),_cat_no);
  }
};


/// \brief average between the leaf obs distros by tree branch length
///
struct prob_gtor_obs_state_time_avg : public prob_gtor_obs {

  typedef prob_gtor_obs_state_time_avg self_t;
  typedef prob_gtor_obs base_t;

  prob_gtor_obs_state_time_avg(const prob_gtor_obs_sml_share& smlm,
                               const unsigned cat_no)
    : base_t(smlm,cat_no) {}

  prob_gtor_obs_state_time_avg(const self_t& s,
                               const prob_gtor_obs_sml_share& smlm)
    : base_t(s,smlm) {}

  virtual self_t* clone(const prob_gtor_obs_sml_share& smlm) const { return new self_t(*this,smlm); }

private:
  virtual void eval_pdistro_update() const {
    _smlm.get_obs_state_time_avg_root(pdistro_mutable(),_cat_no);
  }
};


/// \brief per branch least-squares ('BLS') root
///
struct prob_gtor_obs_root_lsprob : public prob_gtor_obs {

  typedef prob_gtor_obs_root_lsprob self_t;
  typedef prob_gtor_obs base_t;

  prob_gtor_obs_root_lsprob(const prob_gtor_obs_sml_share& smlm,
                            const unsigned cat_no)
    : base_t(smlm,cat_no) {}

  prob_gtor_obs_root_lsprob(const self_t& s,
                            const prob_gtor_obs_sml_share& smlm)
    : base_t(s,smlm) {}

  virtual self_t* clone(const prob_gtor_obs_sml_share& smlm) const { return new self_t(*this,smlm); }

private:
  virtual void eval_pdistro_update() const {
    _smlm.get_root_node_prob(pdistro_mutable(),_cat_no);
  }
};



#ifdef BRANCH_SPECIFIC_OBS
struct prob_gtor_root_lsprob : public prob_gtor_obs {

  typedef prob_gtor_root_lsprob self_t;
  typedef prob_gtor_obs base_t;

  explicit
  prob_gtor_root_lsprob(const prob_gtor_obs_sml_share& smlm)
    : base_t(smlm,0),_node_state_prob(smlm.branch_size()+1,smlm.state_size()) {}

  prob_gtor_root_lsprob(const self_t& s,
                        const prob_gtor_obs_sml_share& smlm)
    : base_t(s,smlm),_node_state_prob(s._node_state_prob) {}

  virtual self_t* clone(const prob_gtor_obs_sml_share& smlm) const { return new self_t(*this,smlm); }

  const prob_t* const * node_state_prob() const {
    return _node_state_prob.val;
  }

private:
  virtual void eval_pdistro_update() const {
    _smlm.get_root_node_prob_nspround(pdistro_mutable(),_cat_no,_node_state_prob.val);
  }

  simple_matrix<prob_t> _node_state_prob;
};
#endif


#endif

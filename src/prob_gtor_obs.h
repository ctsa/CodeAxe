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

// $Id: prob_gtor_obs.h 1122 2008-01-28 23:13:58Z ctsa $

/// \file

#ifndef __PROB_GTOR_OBS_H
#define __PROB_GTOR_OBS_H

#include "observer.h"
#include "prob_gtor.h"
#include "subs_ml_model_comex.h"



/// \brief
///
struct prob_gtor_obs : public prob_gtor, public observer {

  typedef prob_gtor_obs self_t;
  typedef prob_gtor base_t;

  prob_gtor_obs(const prob_gtor_obs_sml_share& smlm,
                const unsigned cat_no)
    : base_t(smlm.state_size()), _is_valid(false), _smlm(smlm), _cat_no(cat_no) {
    observe_notifier(_smlm.get_obs_info_notifier());
  }

  prob_gtor_obs(const self_t& s,
                const prob_gtor_obs_sml_share& smlm)
    : base_t(s), _is_valid(false), _smlm(smlm), _cat_no(s._cat_no) {
    observe_notifier(_smlm.get_obs_info_notifier());
  }

private:
  prob_gtor_obs(const self_t&);

public:

  virtual
  const prob_t* pdistro() const {
    if(! _is_valid) {
      eval_pdistro_update();
      _is_valid=true;
    }
    return base_t::pdistro();
  }

  virtual self_t* clone(const prob_gtor_obs_sml_share&) const = 0;

  // pdistro_update is lazy. calculation is defered until call to pdistro();
  //
  void pdistro_update() { _is_valid=false;}

private:
  virtual void eval_pdistro_update() const = 0;

  virtual
  void recieve_notifier_event(const notifier* n,
                              const EVENT_TYPE::index_t e){
    if(n == &(_smlm.get_obs_info_notifier()) && e == EVENT_TYPE::PDISTRO_CHANGE){
      pdistro_update();
    }
  }


  mutable bool _is_valid;
protected:
  const prob_gtor_obs_sml_share _smlm;
  unsigned _cat_no;
};



#endif

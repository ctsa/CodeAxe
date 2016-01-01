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

// $Id: prob_gtor_param_gc_shifter.h 1121 2008-01-28 22:32:57Z ctsa $

/// \file

#ifndef __PROB_GTOR_PARAM_GC_SHIFTER_H
#define __PROB_GTOR_PARAM_GC_SHIFTER_H


#include "observer.h"
#include "prob_gtor_param.h"
#include "util/bio/bioseq_util.h"

#include <cmath>



inline
smlfloat sigm(const smlfloat x) { return 1./(1.+std::exp(-x)); }

inline
smlfloat inv_sigm(const smlfloat x) { return std::log(x/(1.-x)); }




/// \brief modifies another prob_distro by a single GC-content shifting
/// parameter
///
struct prob_gtor_param_gc_shifter : public prob_gtor_param, public observer {

  typedef prob_gtor_param_gc_shifter self_t;
  typedef prob_gtor_param base_t;

  prob_gtor_param_gc_shifter(const prob_gtor_param& pg,
                             const SITE_MODEL::index_t si)
    : base_t(pg.state_size(),1),_pg(pg), _si(si), _gc_factor(state_size()), _pgdat(state_size()) {

    if(state_size() != SITE_MODEL::state_size(_si)){
      die("inconsistent initialization of _gc_shifter");
    }

    _is_free_block_start[0] = false;
    _is_free_block_stop[0] = false;
    _is_nonnegative[0] = false;

    NUC::index_t nuc[SITE_MODEL::MAX_BASE_SIZE];
    const unsigned ss(state_size());
    const unsigned bsize(SITE_MODEL::base_size(_si));

    for(unsigned i(0);i<ss;++i){
      decode_nuc(_si,i,nuc);
      unsigned gc(0);
      for(unsigned b(0);b<bsize;++b){
        if(nuc[b] == NUC::G || nuc[b] == NUC::C) gc++;
      }
      _gc_factor[i] = (static_cast<smlfloat>(gc)/static_cast<smlfloat>(bsize));
    }

    observe_notifier(_pg);
  }


  prob_gtor_param_gc_shifter(const self_t& s,
                             const prob_gtor_param& pg)
    : base_t(s), _pg(pg), _si(s._si), _gc_factor(s._gc_factor), _pgdat(state_size()) {

    if(state_size() != SITE_MODEL::state_size(_si)){
      die("inconsistent initialization of _gc_shifter");
    }

    observe_notifier(_pg);
  }

  // deal with 1) getting a new copy of _pg and 2) re-registering as an observer of that _pg
  virtual self_t* clone() const { die("unfinished"); return 0; }
  //  virtual self_t* clone(const prob_gtor& pg) const { return new self_t(*this); }

private:
  virtual
  void get_param_from_pdistro() {
    if(! _pgdat.is_init) die("no pgdat init");

    // get gc content of pdistro, psat and psgc, and use this to set param:
    const unsigned ss(state_size());
    prob_t gc_pdistro(0.),gc_psat(0.),gc_psgc(0.);
    for(unsigned i(0);i<ss;++i){
      gc_pdistro+=pdistro()[i]*_gc_factor[i];
      gc_psat+=_pgdat.psat[i]*_gc_factor[i];
      gc_psgc+=_pgdat.psgc[i]*_gc_factor[i];
    }

    const smlfloat target_p((gc_pdistro-gc_psat)/(gc_psgc-gc_psat));
    _param[0]=p_to_param(target_p)-_pgdat.psf;
  }

  void
  get_pgdat_from_pgdistro() {
    /// break starting ps(x) down into two distros based on gc content, so that:
    ///
    /// p(x) = 1/(1-Z)*psat(x)*(1-sig(psf+f))+1/(Z)*psgc(x)*sig(psf+f)
    ///
    /// where f is the shifters free parameter in [-inf,+inf], starting at 0.
    /// get sf,psat and psgc from ps
    ///
    const unsigned ss(state_size());
    const prob_t* pgd(_pg.pdistro());

    smlfloat gcsum(0.);
    for(unsigned i(0);i<ss;++i){
      _pgdat.psat[i]=pgd[i]*(1.-_gc_factor[i]);
      _pgdat.psgc[i]=pgd[i]*_gc_factor[i];
      gcsum += _pgdat.psgc[i];
    }

    for(unsigned i(0);i<ss;++i){
      _pgdat.psat[i] /= (1.-gcsum);
      _pgdat.psgc[i] /= gcsum;
    }

    _pgdat.psf = p_to_param(gcsum);
    _pgdat.is_init=true;
  }

  virtual
  void get_pdistro_from_param() {
    if(! _pgdat.is_init) die("no pgdat init");
    const unsigned ss(state_size());

    const prob_t pgc(param_to_p(_pgdat.psf+_param[0]));

    for(unsigned i(0);i<ss;++i) {
      pdistro()[i] = _pgdat.psat[i]*(1.-pgc)+_pgdat.psgc[i]*pgc;
    }
  }

  virtual
  void param_norm() {}

  virtual
  void recieve_notifier_event(const notifier* n,
                              const EVENT_TYPE::index_t e){
    if(n == &_pg && e == EVENT_TYPE::PDISTRO_CHANGE){
      get_pgdat_from_pgdistro();
    }
  }

  static
  smlfloat p_to_param(prob_t x) { return inv_sigm(x); }

  static
  prob_t param_to_p(smlfloat x) { return sigm(x); }


  struct pgdat_t {

    pgdat_t(const unsigned s) : psat(s), psgc(s), psf(0.), is_init(false) {}

    std::vector<prob_t> psat;
    std::vector<prob_t> psgc;
    smlfloat psf;
    bool is_init;
  };

  const prob_gtor_param& _pg;
  const SITE_MODEL::index_t _si;
  std::vector<smlfloat> _gc_factor;
  pgdat_t _pgdat;
};


#endif

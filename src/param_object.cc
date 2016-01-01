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

// $Id: param_object.cc 1113 2008-01-27 00:23:59Z ctsa $

/// \file

#include "param_object.h"
#include "param_util.h"
#include "util/general/die.h"



unsigned
param_object::
param_size(const PARAM_VIEW::index_t pv) const {
  const unsigned ps(_param.size());
  if       (pv == PARAM_VIEW::ALL){
    return ps;
  } else if(pv == PARAM_VIEW::MIN_PARAM ||
            pv == PARAM_VIEW::TRAINABLE){
    unsigned pst(0);
    for(unsigned i(0);i<ps;++i) if(_is_train_param[i]) pst++;
    if(pv == PARAM_VIEW::TRAINABLE) { return pst; }
    else                            { return pst+pseudo_param_size(); }

  } else if(pv == PARAM_VIEW::INDY_MIN_PARAM){
    unsigned pst(0);
    bool is_in_block(false);
    bool is_first_in_block(false);
    for(unsigned i(0);i<ps;++i) {
      if(_is_free_block_start[i]) {
        is_in_block=true;
        is_first_in_block=true;
      }

      if(_is_train_param[i]){
        if(is_in_block && is_first_in_block){
          is_first_in_block=false;
        } else { pst++; }
      }

      if(_is_free_block_stop[i]) is_in_block=false;
    }
    return pst;
  } else { die("param_size: unknown param view"); }
}



smlfloat*
param_object::
param_state(smlfloat* x,
            const PARAM_VIEW::index_t pv) const {
  const unsigned ps(param_size());

  if       (pv == PARAM_VIEW::ALL){
    for(unsigned i(0);i<ps;++i) {*x =_param[i]; ++x;}

  } else if(pv == PARAM_VIEW::TRAINABLE){
    for(unsigned i(0);i<ps;++i){
      if(_is_train_param[i]) {*x=_param[i]; ++x;}
    }

  } else if(pv == PARAM_VIEW::MIN_PARAM){
    const unsigned pps(pseudo_param_size());
    for(unsigned i(0);i<pps;++i) { *x = 1.; ++x; }

    for(unsigned i(0);i<ps;++i){
      if(_is_train_param[i]) {
        *x = trainable_param_transform(i);
        ++x;
      }
    }

  } else if(pv == PARAM_VIEW::INDY_MIN_PARAM){
    bool is_in_block(false);
    for(unsigned i(0);i<ps;++i){
      if(_is_free_block_start[i]) is_in_block=true;

      if(_is_train_param[i]) {
        if(is_in_block){
          if(! _is_nonnegative[i]) die("invalid block flag combo in param_object");
        } else {
          *x = trainable_param_transform(i);
          ++x;
        }
      }

      if(_is_free_block_stop[i]) {
        is_in_block=false;

        x=dependent_to_independent_param_block(_is_train_param,
                                               _is_free_block_start,
                                               i,
                                               _param,
                                               x);
      }
    }

  } else { die("param_state: unknown param view"); }

  return x;
}



const smlfloat*
param_object::
set_param_state(const smlfloat* x,
                const PARAM_VIEW::index_t pv) {

  const unsigned ps(param_size());

  if       (pv == PARAM_VIEW::ALL){
    for(unsigned i(0);i<ps;++i) {
      _param[i]=*x;
      if(_is_nonnegative[i] && _param[i] < 0.) die("invalid param value");
      ++x;
    }

  } else if(pv == PARAM_VIEW::TRAINABLE){
    for(unsigned i(0);i<ps;++i) {
      if(_is_train_param[i]) {
        _param[i] = *x;
        if(_is_nonnegative[i] && _param[i] < 0.) die("invalid param value");
        ++x;
      }
    }

  } else if(pv == PARAM_VIEW::MIN_PARAM){
    //    typedef typename std::iterator_traits<InputIterator>::value_type InputT;
    typedef smlfloat InputT;

    const unsigned pps(pseudo_param_size());
    std::vector<InputT> p(pps);

    // 1) store pseudoparams
    for(unsigned i(0);i<pps;++i) {
      p[i] = *x;
      ++x;
    }

    // 2) read regular params
    for(unsigned i(0);i<ps;++i) {
      if(_is_train_param[i]) {
        trainable_param_untransform(i,(_is_nonnegative[i] ? std::fabs(*x) : *x ));
        ++x;
      }
    }

    // 3) go back to pseudoparams after regular params have been read:
    for(unsigned i(0);i<pps;++i) { set_pseudo_param(i,p[i]); }

  } else if(pv == PARAM_VIEW::INDY_MIN_PARAM){

    x=set_indy_min_param(_is_train_param,
                         _is_free_block_start,
                         _is_free_block_stop,
                         _is_nonnegative,
                         ps,
                         _param,
                         x);

    for(unsigned i(0);i<ps;++i){
      if(_is_train_param[i]) trainable_param_untransform(i,_param[i]);
    }

  } else { die("set_param_state: unknown param view"); }

  post_write_param_state();
  return x;
}



void
param_object::
set_param_state(const smlfloat x) {

  const unsigned ps(param_size());
  for(unsigned i(0);i<ps;++i) {
    _param[i]=x;
    if(_is_nonnegative[i] && _param[i] < 0.) die("invalid param value");
  }

  post_write_param_state();
}

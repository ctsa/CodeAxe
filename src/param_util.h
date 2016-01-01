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

// $Id: param_util.h 791 2007-09-06 21:55:22Z ctsa $

/// \file
///
/// \brief functions to convert parameter sets between dependent and
/// independent forms
///
/// Dependent sets of parameters, such as N parameters used to express
/// the probability distribution of N states, are used internally in
/// the model (1) for ease of interpretation and (2) as input to the
/// conjugate direction minimizer. These must be converted to
/// independent sets for any minimizer which takes the derivative of
/// the score wrt the parameter vector (conj grad, bfgs...). For best
/// minimization, it is better if we can give the minimizer values
/// which can be set over a large range (nonnegative) -- this makes
/// passing all parameters but one of a prob distribution undesirable,
/// b/c hacky things need to be done to keep the distro sum == 1
///
/// A further complication for the sml model is that any parameter
/// within the dependent set may be locked, and thus a subset of the
/// dependent parameter set needs to be transformed back and forth to
/// the independent set.
///
/// A naive solution is to express the unlocked parameters in the
/// dependent set as a ratio of each value to, for instance, the first
/// unlocked value in the distribution. However this restricts the
/// denominator parameter to be nonzero.
///
/// The following functions implement a more robust transformation, assuming
/// that all dependent parameter values are nonnegative and the sum of all
/// trainable dependent parameters is nonzero:
///
/// Definitions:
///
/// \f$D = \{ d_x \}\f$ is the set of \f$X\f$ dependent params, \f$d_x\f$ (excluding locked params)
///
/// \f$I = \{ i_x \}\f$ is the set of \f$X-1\f$ independent params, \f$i_x\f$
///
///
/// D to I:
///
///   \f$a = {sum_{x}{d_x}} \over X\f$
///
///   \f$bal = (d_1/a)+1\f$
///
///   \f$i_{x-1} = ((d_x/a)+1)/bal , x=2,X\f$
///
///
/// I to D:
///
///   \f$a = 2X \over {1+\sum_{x=2,X}{i_{x-1}}}\f$
///
///   \f$d_1 = a-1\f$
///
///   \f$d_x = (i_{x-1}a)-1 , x=2,X\f$
///


#ifndef __PARAM_UTIL_H
#define __PARAM_UTIL_H

#include "subs_ml_types.h"
#include "util/general/die.h"

#include <cassert>
#include <cmath>

#include <vector>

template <typename OutputIterator>
OutputIterator
dependent_to_independent_param_block(const std::vector<bool>& is_train_param,
                                     const std::vector<bool>& is_free_block_start,
                                     const unsigned i,
                                     const std::vector<smlfloat>& param,
                                     OutputIterator x){

  smlfloat block_scale(0);
  unsigned block_size(0);
  unsigned block_start(0);

  for(unsigned j(i);;--j){
    if(is_train_param[j]){
      block_scale += param[j];
      block_size++;
    }
    if(is_free_block_start[j]) {
      block_start=j;
      break;
    }
    if(j==0) die("Invalid independent param block structure.");
  }

  if(block_size==0) return x;

  const smlfloat A(block_scale/static_cast<smlfloat>(block_size));

  smlfloat balance_param(0);

  bool is_first(true);
  for(unsigned j(block_start);j<=i;++j){
    if(is_train_param[j]){
      if(is_first){
        balance_param = ((param[j]/A)+1);
        is_first = false;
      } else {
        *x = ((param[j]/A)+1)/balance_param;
        ++x;
      }
    }
  }
  return x;
}


void
independent_to_dependent_param_block(const std::vector<bool>& is_train_param,
                                     const std::vector<bool>& is_free_block_start,
                                     const unsigned i,
                                     std::vector<smlfloat>& param);


/// \brief does the indy->dep param block transformation over the
/// whole set of parameters given to the minimizer
///
template <typename InputIterator>
InputIterator
set_indy_min_param(const std::vector<bool>& is_train_param,
                   const std::vector<bool>& is_free_block_start,
                   const std::vector<bool>& is_free_block_stop,
                   const std::vector<bool>& is_nonnegative,
                   const unsigned ps,
                   std::vector<smlfloat>& param,
                   InputIterator x){

  bool is_in_block(false);
  bool is_first_in_block(false);

  for(unsigned i(0);i<ps;++i){
    if(is_free_block_start[i]){
      is_in_block=true;
      is_first_in_block=true;
    }

    if(is_train_param[i]){
      if(is_in_block && is_first_in_block){
        is_first_in_block=false;
        param[i] = 1.;
      } else {
        param[i] = *x;
        if(is_nonnegative[i]) param[i] = std::fabs(param[i]);
        ++x;
      }
    }

    if(is_free_block_stop[i]) {
      assert(is_in_block);
      is_in_block=false;

      independent_to_dependent_param_block(is_train_param,
                                           is_free_block_start,
                                           i,
                                           param);
    }
  }
  return x;
}


#endif

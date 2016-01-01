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

// $Id: random_util.cc 743 2007-08-14 15:47:12Z ctsa $

/// \file
///


#include "random_util.h"

#include <ctime>


#ifdef DEBUG
#include "../general/log.h"
#include <ostream>
#endif

#ifdef USE_BOOST_RANDOM
// using static rng here is suboptimal, but it replicates the behavior of the
// c standard library random functions:
base_generator_type br_rangen(0u);
br_uni_func_type br_uni_func(br_rangen,boost::uniform_real<>());
#endif

long int
random_init(long int seed){
  if(seed==-1) seed = time(0);
#ifdef USE_BOOST_RANDOM
  seed=std::abs(seed);
  br_rangen.seed(static_cast<unsigned>(seed));
#else
  srand48(seed);
#endif
#ifdef DEBUG
  log_os << "Initialized rng: seed= " << seed << "\n";
#endif
  return seed;
};

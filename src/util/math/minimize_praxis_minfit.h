// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// CodeAxe : phylogenetic analysis and simulation tools
//
//   http://www.phrap.org
//
//

// $Id: minimize_praxis_minfit.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __MINIMIZE_PRAXIS_MINFIT_C
#define __MINIMIZE_PRAXIS_MINFIT_C

namespace PRAXIS {

/* singular value decomposition */
template <typename FloatType>
void
minfit(const int n,
       FloatType eps,
       const FloatType tol,
       FloatType* ab,
       FloatType* q,
       FloatType* workspace);

}

#include "minimize_praxis_minfit.hh"

#endif

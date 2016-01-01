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

// $Id: linear_solver.cc 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifdef USE_LAPACK

#include "clapack_util.h"
#include "linear_solver.h"
#include "math_util_exception.h"



void
linear_solver(const double* A,
              const double* b,
              double* x,
              const unsigned N){

  using namespace CLAPACK;

  integer one(1);
  integer n(static_cast<integer>(N));
  integer info(0);

  integer* iworkspace(new integer[N]);
  integer* pivots(iworkspace);

  double* dworkspace(new double[N*N]);
  double* Acopy(dworkspace);

  std::copy(b,b+N,x);
  std::copy(A,A+N*N,Acopy);

  dgesv_(&n,&one,Acopy,&n,pivots,x,&n,&info);

  if(info){ throw math_util_exception("gesvx (LAPACK linear solver) failure"); }

  delete[] iworkspace;
  delete[] dworkspace;
}



void
linear_solver_refined(const double* A,
                      const double* b,
                      double* x,
                      double& ferr,
                      const unsigned N){

  using namespace CLAPACK;

  integer one(1);
  integer n(static_cast<integer>(N));
  integer info(0);
  double rcond(0.);
  double berr(0.);
  char equed[2] = "X";

  integer* iworkspace(new integer[2*N]);
  integer* iwork(iworkspace);
  integer* pivots(iwork+N);

  double* dworkspace(new double[(7+2*N)*N]);
  double* work(dworkspace);
  double* row_equid(work+4*N);
  double* col_equid(row_equid+N);
  double* bcopy(col_equid+N);
  double* Acopy(bcopy+N);
  double* Acopy2(Acopy+N*N);

  std::copy(b,b+N,bcopy);
  std::copy(A,A+N*N,Acopy);

  dgesvx_("E","N",&n,&one,Acopy,&n,Acopy2,&n,pivots,equed,row_equid,col_equid,
          bcopy,&n,x,&n,&rcond,&ferr,&berr,work,iwork,&info);

  if(info){ throw math_util_exception("gesvx (LAPACK linear solver) failure"); }

  delete[] iworkspace;
  delete[] dworkspace;
}

#endif

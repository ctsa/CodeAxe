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

// $Id: stationary_pdistro.cc 1057 2007-12-10 23:31:52Z ctsa $

/// \file

#include "simple_util.h"
#include "stationary_pdistro.h"
#include "subs_ml_ptol.h"
#include "substk_exception.h"
#include "util/general/die.h"
#include "util/general/log.h"
#include "util/math/clapack_util.h"
#include "util/math/matrix_util.h"
#include "util/math/matrix_util_io.h"
#include "util/math/prob_util.h"
#include "util/math/prob_util_io.h"

#include <cmath>

#include <algorithm>
#include <iomanip>
#include <ostream>



const smlfloat SD_PTOL(1.e-5);

const smlfloat STAT_ZERO_THRESH(1.e-11);


#ifdef USE_LAPACK
/// LAPACK LU decomposition:
///
static
void
get_matrix_LU(const unsigned N,
              double* mat,
              int* pivots = 0){

  using namespace CLAPACK;

  integer n(static_cast<integer>(N));

  integer* local_pivots(new integer[N]);

  integer info(0);

  dgetrf_(&n,&n,mat,&n,local_pivots,&info);

  if( info < 0 ){ die("getrf (gen lu decomp) failure."); }
  // if info is *greater* than 0, then the matrix is singular,
  // allowed in this case...

  if(pivots) std::copy(local_pivots,local_pivots+N,pivots);

  delete [] local_pivots;
}
#endif



void
get_stationary_pdistro_from_rates(smlfloat* stat_pdistro,
                                  const smlfloat* rates,
                                  const unsigned N,
                                  const bool
#ifndef USE_LAPACK
                                  is_start_from_input
#endif
                                  ){

#ifdef USE_LAPACK

#ifdef DEBUG
  for(unsigned i(0);i<N;++i){
    smlfloat sum(0.);
    for(unsigned j(0);j<N;++j){
      if(i!=j) {
        smlfloat val(rates[j+i*N]);
        assert(val>=0);
        sum -= val;
      }
    }
    assert(std::fabs(sum-rates[i+i*N])<STAT_ZERO_THRESH);
  }
#endif

  // solution as discussed in page 3 of Stewart, WJ "Numerical methods
  // for computing stationary distribution of finite irreducible Markov
  // chains"
  //
  double* matrix(new double[N*N]);

  std::copy(rates,rates+N*N,matrix);

  get_matrix_LU(N,matrix);

  // arbitrarily set last value to 1. (any constant will do)
  stat_pdistro[N-1] = 1.;

  for(int i(N-2);i>=0;--i){
    stat_pdistro[i] = 0.;
    for(int j(N-1);j>i;--j){
      stat_pdistro[i] -= stat_pdistro[j]*matrix[i+j*N];
    }
    stat_pdistro[i] /= matrix[i+i*N];
  }
  delete [] matrix;

  pdistro_norm(stat_pdistro,stat_pdistro+N);
#ifdef DEBUG
  pdistro_check(stat_pdistro,N,SUBS_ML_PTOL);

  smlfloat* testvec(new smlfloat[N]);
  array_zero(testvec,N);
  matrix_vector_mult_sum(rates,stat_pdistro,testvec,N,N);
  for(unsigned i(0);i<N;++i){
    if(std::fabs(testvec[i])>STAT_ZERO_THRESH){
      log_os << "pdistro_test failed! i,val " << i << " " << testvec[i] << "\n";
    }
  }

  delete [] testvec;
#endif
#else
  /// \todo write a non-lapack fallback lu decomposer, then get rid of this code:
  static const smlfloat MIN_TIME_DELTA(2.e-2);

  smlfloat norm(0.);
  for(unsigned i(0);i<N;++i){
    for(unsigned j(0);j<N;++j) if(i!=j) norm += rates[i+j*N];
  }

  norm /= static_cast<smlfloat>(N);

  const smlfloat rate_factor(MIN_TIME_DELTA/norm);

  smlfloat* pstep = new smlfloat[N*N];
  for(unsigned i(0);i<N;++i){
    smlfloat sum(0.);
    for(unsigned j(0);j<N;++j){
      if(i==j) continue;
      pstep[j+i*N] = rates[j+i*N]*rate_factor;
      sum += pstep[j+i*N];
    }
    pstep[i+i*N] = 1.-sum;
  }

  get_stationary_pdistro_from_pmatrix(stat_pdistro,pstep,N,is_start_from_input);

  delete [] pstep;
#endif
}


void
get_stationary_pdistro_from_pmatrix(smlfloat* stat_pdistro,
                                    const smlfloat* pstep,
                                    const unsigned N,
                                    const bool is_start_from_input){

  /// \todo replace iterative solution with principle eigenvector of pstep:
  ///
  static const unsigned MAX_ITER(1000000);
  static const smlfloat STAT_CONVERGE_THRESH(1.e-8);
  static const unsigned MIN_FREEZE(10);

  static const unsigned max_peeksize(4);
  const unsigned peeksize(std::min(N,max_peeksize));

#ifdef DEBUG
  for(unsigned i(0);i<N;++i){
    pdistro_check(pstep+(i*N),N,SD_PTOL);
  }
#endif

  // setup starting distro:
  if(! is_start_from_input) pdistro_unif(stat_pdistro,stat_pdistro+N);

  simple_array<smlfloat> stat_pdistro2_a(N);
  smlfloat* stat_pdistro2(stat_pdistro2_a.ptr());

  // stopping data:
  unsigned total_freeze(0);
  unsigned last_freeze(0);

  unsigned iter(0);
  for(;iter<MAX_ITER;iter++){
    array_zero(stat_pdistro2,N);
    matrix_vector_mult_sum(pstep,stat_pdistro,stat_pdistro2,N,N);

    if((iter+10000)>MAX_ITER){
      log_os << std::setprecision(12);
      for(unsigned i(0);i<peeksize;++i) {
        if(i) log_os << " ";
        log_os << stat_pdistro2[i];
      }
      log_os << "\n";
    }

    // test for finish:
    unsigned j(0);
    for(;j<N;++j) if(std::fabs(stat_pdistro[j]-stat_pdistro2[j])>STAT_CONVERGE_THRESH) break;
    if(j==N) {
      if(last_freeze==(iter-1)) total_freeze++;
      else                      total_freeze=1;
      if(total_freeze>=MIN_FREEZE) break;
      last_freeze=iter;
    }
    std::swap(stat_pdistro,stat_pdistro2);
  }

  if(iter==MAX_ITER){
    std::ostringstream oss;
    oss << "stat_distro_from_pmatrix(): stationary distribution calc failed to converge\n";
    for(unsigned i(0);i<peeksize;++i){
      for(unsigned j(0);j<peeksize;++j){
        oss << pstep[j+i*N] << " ";
      }
      oss << "\n";
    }
    throw substk_exception(oss.str().c_str());
  }
#if DEBUG
  pdistro_check(stat_pdistro,N,SD_PTOL);
  log_os << "Total stat iterations: " << iter << "\n";
#endif
}


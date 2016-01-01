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

// $Id: subs_ml_model_root_from_obs.cc 1160 2008-03-20 19:18:02Z ctsa $

/// \file

#include "bi_tree.h"
#include "cat_manager.h"
#include "rate_gtor.h"
#include "time_gtor.h"
#include "simple_util.h"
#include "subs_ml_model_min_options.h"
#include "subs_ml_model_root_from_obs.h"
#include "subs_ml_ptol.h"
#include "subs_prob_util.h"
#include "util/general/log.h"
#include "util/math/array_util.h"
#include "util/math/linear_solver.h"
#include "util/math/matrix_util.h"

#ifdef BRANCH_SPECIFIC_OBS
#include "rate_edge_dependencies.h"
#endif

// conj gradient minimization of the least-sq root problem
// is compiled in by default, alternatives are NNLS and conj
// direction minimizer
//
#ifdef USE_NNLS_ROOT
#include "util/math/nnls.h"
#else

#ifndef USE_CONJ_DIR_ROOT
#include "util/math/minimize_conj_gradient.h"
#else
#include "util/math/minimize_conj_direction.h"
#endif

#endif

#include "util/math/minfunc_interface.h"
#include "util/math/prob_util.h"

#include <iomanip>
#include <ostream>


//#define DEBUG_BGR


static
void
tree_branch_down_obs_state_time_avg(const bi_tree_node* b,
                                    const prob_t*& branch_distro,
                                    prob_t** down_distro_tmp,
                                    const time_gtor& tgm,
                                    const unsigned n_states,
                                    const prob_t* const * leaf_distro){

  const unsigned n_time_cats(tgm.time_cat_size());

  simple_array<prob_t> time_cat_pdistro(n_time_cats);
  tgm.time_cat_pdistro(time_cat_pdistro.ptr());


  if(b==0){
    die("tree_branch_down_obs_state_time_avg: null tree node");

  } else if(b->is_leaf()){  // leaf node
    const unsigned leaf_id(b->leaf_id());
    branch_distro = leaf_distro[leaf_id];

  } else {  // assume child2 exists

    const unsigned c1i(b->child1()->branch_id());
    const unsigned c2i(b->child2()->branch_id());

    smlfloat c1ntime(0.);
    smlfloat c2ntime(0.);

    for(unsigned tc(0);tc<n_time_cats;++tc){
      const smlfloat c1time(tgm.time_cat_branch_time(c1i,tc));
      const smlfloat c2time(tgm.time_cat_branch_time(c2i,tc));
      const smlfloat c1ntime_cat(c1time/(c1time+c2time));
      const smlfloat c2ntime_cat(c2time/(c1time+c2time));

      const prob_t time_cat_prob(time_cat_pdistro[tc]);
      c1ntime += time_cat_prob*c1ntime_cat;
      c2ntime += time_cat_prob*c2ntime_cat;
    }

    const prob_t* child1_distro(0);
    const prob_t* child2_distro(0);

    tree_branch_down_obs_state_time_avg(b->child1(),child1_distro,down_distro_tmp,
                                        tgm,n_states,leaf_distro);
    tree_branch_down_obs_state_time_avg(b->child2(),child2_distro,down_distro_tmp,
                                        tgm,n_states,leaf_distro);

    const unsigned node_id(b->node_id());
    prob_t* branch_distro_tmp(down_distro_tmp[node_id]);
    branch_distro = branch_distro_tmp;

    for(unsigned i(0);i<n_states;++i){
      branch_distro_tmp[i] += child1_distro[i]*c1ntime + child2_distro[i]*c2ntime;
    }
  }
}



void
subs_ml_model_obs_state_time_avg_root(const subs_ml_model& mdl,
                                      prob_t* root_down_prob,
                                      const unsigned cat_no){

  const time_gtor& tgm(mdl.get_time_gtor());
  const bi_tree& tree(mdl.tree());
  const unsigned n_states(mdl.state_size());

  const prob_t* const * leaf_distro(mdl.obs_seq_state_distro_cat(cat_no));

  simple_matrix<prob_t> down_distro_tmp(tree.node_size(),n_states);

  const prob_t* root_distro(0);
  tree_branch_down_obs_state_time_avg(tree.root(),root_distro,down_distro_tmp.ptr(),
                                      tgm,n_states,leaf_distro);

  for(unsigned i(0);i<n_states;++i){ root_down_prob[i] = root_distro[i]; }
}






struct parent_prob_leastsq_minfunc_base : public minfunc_gradient_interface<smlfloat> {

  typedef parent_prob_leastsq_minfunc_base self_t;

  parent_prob_leastsq_minfunc_base(const prob_t * const tprob,
                                   const prob_t * const child_distro,
                                   const unsigned n_states)
    : _tprob(tprob), _child_distro(child_distro), _n_states(n_states),
      _parent_prob(_n_states), _tmp_store(_n_states) {}

private:
  parent_prob_leastsq_minfunc_base(const self_t&);
  self_t& operator=(const self_t&);

public:

  virtual
  void p_to_parent_prob(const smlfloat* p,
                        prob_t* parent_prob) const = 0;

  virtual
  void parent_prob_to_p(const prob_t* parent_prob,
                        smlfloat* p) = 0;

  virtual smlfloat val(const smlfloat* p){
    return val_int(p);
  }

  // i=parent_state
  // j=leaf_state
  // x=parent_param
  // t1(x,j) = sum_i( x(i)* tprob(j,i) )
  // t2(x,j) = t1(x,j) - ldistro(j)
  // cost(x)= sum_j( t2(x,j)**2 )
  //
  smlfloat val_int(const smlfloat* p,
                   smlfloat* tval = 0) {

    p_to_parent_prob(p,_parent_prob.ptr());

    smlfloat cost(0.);
    smlfloat* tval_int(0);

    if(tval) tval_int=tval;
    else     tval_int=_tmp_store.ptr();

    for(unsigned j(0);j<_n_states;++j){
      tval_int[j] = array_dot(_parent_prob.ptr(),_tprob+j*_n_states,_n_states);
    }

    for(unsigned j(0);j<_n_states;++j){
      tval_int[j] -=_child_distro[j];
      cost += tval_int[j]*tval_int[j];
    }

    cost *= cost_scale_factor();
    return cost;
  }

  // double sided numerical approx derivative:
  virtual
  smlfloat dval(const smlfloat* p,
                smlfloat* dp){

    static const smlfloat delta(1e-5);
    static const smlfloat deltamin(1e-7);

    const unsigned ndim(dim());

    const smlfloat start_cost(val(p));

    simple_array<smlfloat> vcopy(ndim);
    std::copy(p,p+ndim,vcopy.ptr());
    for(unsigned i(0);i<ndim;++i){
      const smlfloat delchunk(std::max(std::fabs(vcopy[i]*delta),deltamin));
      vcopy[i] = p[i]+delchunk;
      const smlfloat pluscost(val(vcopy.ptr()));
      vcopy[i] = p[i]-delchunk;
      const smlfloat minuscost(val(vcopy.ptr()));
      dp[i] = (pluscost-minuscost)/(delchunk*2.);
    }

    return start_cost;
  }

protected:
  smlfloat cost_scale_factor() const {
    static const smlfloat s(1e0);
    return s/(_n_states);
  }

  const prob_t * const _tprob;
  const prob_t * const _child_distro;
  const unsigned _n_states;
  simple_array<prob_t> _parent_prob;
  simple_array<prob_t> _tmp_store;
};



struct parent_prob_leastsq_minfunc : public parent_prob_leastsq_minfunc_base {

  typedef parent_prob_leastsq_minfunc_base base_t;

  parent_prob_leastsq_minfunc(const prob_t * const tprob,
                              const prob_t * const child_distro,
                              const unsigned n_states)
    : base_t(tprob,child_distro,n_states) {}

  virtual
  unsigned dim() const { return _n_states; }

  virtual
  void
  p_to_parent_prob(const smlfloat* p,
                   prob_t* parent_prob) const {

    const unsigned ndim(dim());

    for(unsigned i(0);i<ndim;++i){
      parent_prob[i] = std::fabs(p[i]);
    }
    pdistro_norm(parent_prob,parent_prob+_n_states);
  }

  virtual
  void
  parent_prob_to_p(const prob_t* parent_prob,
                   smlfloat* p) {

    const unsigned ndim(dim());

    for(unsigned i(0);i<ndim;++i){ p[i] = parent_prob[i]; }
  }
};



struct parent_prob_leastsq_minfunc_indy : public parent_prob_leastsq_minfunc_base {

  typedef parent_prob_leastsq_minfunc_base base_t;

  parent_prob_leastsq_minfunc_indy(const prob_t * const tprob,
                                   const prob_t * const child_distro,
                                   const unsigned n_states)
    : base_t(tprob,child_distro,n_states), _parent_key(0) {}

  virtual
  unsigned dim() const { return _n_states-1; }

  virtual
  void
  p_to_parent_prob(const smlfloat* p,
                   prob_t* parent_prob) const {

    const unsigned ndim(dim());

    for(unsigned i(0);i<ndim;++i){
      parent_prob[i] = std::fabs(p[i])*(_parent_key+DIV_ZERO_PROTECT_OFFSET());
    }
    parent_prob[ndim] = _parent_key;
    pdistro_norm(parent_prob,parent_prob+_n_states);
  }


  virtual
  void
  parent_prob_to_p(const prob_t* parent_prob,
                   smlfloat* p) {

    const unsigned ndim(dim());

    _parent_key=parent_prob[ndim];

    for(unsigned i(0);i<ndim;++i){
      p[i] = parent_prob[i]/(_parent_key+DIV_ZERO_PROTECT_OFFSET());
    }
  }


  virtual
  smlfloat dval(const smlfloat* p,
                smlfloat* dp){

    // derivation of analytical derivative:
    //
    // key:
    //
    // i=parent_state
    // j=leaf_state
    // k-model_param_index
    //
    // x=parent_prob
    // y=model_param
    // pk=model_param_transform_key
    // DZ=zero_protect_offset
    // w=pk+DZ
    //
    //
    // parent_prob_function (assuming all y(k)>=0.):
    //
    // x(i!=N,y) = y(i)*w/(sum_k(p(k))*w+pk )
    // x(i!=N,y) = y(i)/(sum_k!=i(y(k))+y(i)+pk/w)
    // AA=sum_k(y(k))+pk/w
    // x(i!=N,y) = y(i)/(AA)
    //
    // x(i==N,y) = pk/(sum_k( y(k)*w )+pk )
    // x(i==N,y) = (pk/w)/(AA)
    //
    //
    // parent_prob_function:
    //
    // x(i!=N,y(k==i)>=0) =  y(i)/(sum_k!=i(std::abs(y(k))+y(i)+pk/w )
    // x(i!=N,y(k==i)<0)  = -y(i)/(sum_k!=i(std::abs(y(k))-y(i)+pk/w )
    //
    // x(i!=N,y(k==z,z!=i)>=0) = std::abs(y(i))/(sum_k!=z(std::abs(y(k))+y(z)+pk/w )
    // x(i!=N,y(k==z,z!=i)<0)  = std::abs(y(i))/(sum_k!=z(std::abs(y(k))-y(z)+pk/w )
    //
    // x(i==N,y(k==z,z!=i)>=0) = (pk/w)/(sum_k!=z(y(k))+y(z)+pk/w )
    // x(i==N,y(k==z,z!=i)<0)  = (pk/w)/(sum_k!=z(y(k))-y(z)+pk/w )
    //
    //
    // parent_prob_derivative (assuming all y(k)>=0.):
    //
    // dx(i,y)_dy(k) =
    //   dx(i!=N,y)_dy(k==i) = 1/(AA) - y(i)/(AA)**2
    //   dx(i!=N,y)_dy(k!=i) = - y(i)/(AA)**2
    //   dx(i==N,y)_dy(k!=i) = - (pk/w)/(AA)**2
    //
    //
    // parent_prob_derivative:
    //
    // dx(i,y)_dy(k),y(k)<0 =
    //   dx(i!=N,y)_dy(k==i) = - 1/(AA) - y(i)/(AA)**2  !!!note the exception here!!!
    //   dx(i!=N,y)_dy(k!=i) = y(i)/(AA)**2
    //   dx(i==N,y)_dy(k!=i) = (pk/w)/(AA)**2
    //
    //
    // cost_function:
    //
    // t1(x,j) = sum_i( x(i)* tprob(j,i) )
    // t2(x,j) = t1(x,j) - ldistro(j)
    // cost(x) = sum_j( t2(x,j)**2 )
    //
    //
    // cost_derivative:
    //
    // dt1(x,j)_dx(i) = tprob(j,i)*dx(i)
    // dt2(x,j)_dx(i) = tprob(j,i)*dx(i)
    // dcost(x)_dx(i) = sum_j( 2*t2(x,j)*tprob(j,i)*dx(i) )
    // dcost(x)_dx(i) = sum_j( t2(x,j)*tprob(j,i) )*2.*dx(i)
    //
    // dcost(y)_dy(k) = sum_i( dcost(x)_dx(i)*dx(i,y)_dy(k) )
    //

    const unsigned ndim(dim());

    simple_array<smlfloat> tval(_n_states);

    const smlfloat start_cost(val_int(p,tval.ptr()));

    smlfloat psum(0.);
    for(unsigned i(0);i<ndim;++i) psum += std::fabs(p[i]);

    const smlfloat Aconst(psum+1/(1.+DIV_ZERO_PROTECT_OFFSET()/_parent_key));
    const smlfloat inv_Aconst(1./Aconst);
    const smlfloat inv_Aconst2(inv_Aconst*inv_Aconst);

    const smlfloat sum_scale(2.*cost_scale_factor());

    std::vector<smlfloat> dcost_dparent(_n_states,0.);

    for(unsigned i(0);i<_n_states;++i){
      for(unsigned j(0);j<_n_states;++j){
        dcost_dparent[i] += tval[j]*_tprob[j*_n_states+i];
      }
    }

    for(unsigned i(0);i<ndim;++i){
      const smlfloat absp(std::fabs(p[i]));

      dp[i] = 0.;
      for(unsigned ii(0);ii<_n_states;++ii){
        smlfloat dparent_ii_dp;
        if       (ii==i) {
          dparent_ii_dp = (inv_Aconst)*(1.-p[i]*inv_Aconst); // intentionally *not* abs(p)
        } else if(ii==ndim) {
          dparent_ii_dp = -_parent_key*inv_Aconst2/(_parent_key+DIV_ZERO_PROTECT_OFFSET());
        } else {
          dparent_ii_dp = -std::abs(p[ii])*inv_Aconst2;
        }
        dp[i] += dcost_dparent[ii]*dparent_ii_dp;
      }

      dp[i] *= sum_scale;
      if(p[i]<0.) { dp[i] *= -1; }

      // smooth deriv as abs(p[i]) approaches 0, but only if deriv is descending towards 0.:
      const smlfloat absdp(std::fabs(dp[i]));
      if(absdp>0.){
        const smlfloat edge_dist(1./(_n_states));
        const smlfloat gg( _parent_prob[i]/edge_dist);
        const smlfloat ff(absp/(10.*absdp));
        //      log_os << "i,p,prob,dp,gg,ff " << i << " " << p[i] << " " << _parent_prob.val[i] << " " << dp[i] << " " << gg << " " << ff << "\n";
        if(gg < 1. && ff < 1. && (dp[i]*p[i] > 0. || (absp <= 0. && dp[i] > 0.))){
          //log_os << "edge,i " << i << "\n";
          dp[i] *= gg;
        }
        //      log_os << "dval2 i,dp: " << i << " " << dp[i] << "\n";
      }
    }

    return start_cost;
  }


private:
  static double DIV_ZERO_PROTECT_OFFSET(){
    static const double d(0.1);
    return d;
  }

  smlfloat _parent_key;
};



struct tbdp_info {
  tbdp_info(bool is_invalid_prob_init,
            prob_t* tmp_store_init,
            prob_t** const branch_down_prob_storage_init,
            const time_gtor& tgm_init,
            const unsigned n_states_init,
            const prob_t* const * const leaf_distro_init,
            const prob_t* const * const branch_tprob_init,
            prob_t** node_state_prob_init)
    : is_invalid_prob(is_invalid_prob_init),
      tmp_store(tmp_store_init),
      branch_down_prob_storage(branch_down_prob_storage_init),
      tgm(tgm_init),
      n_states(n_states_init),
      leaf_distro(leaf_distro_init),
      branch_tprob(branch_tprob_init),
      tprob_tmp(n_states*n_states),
      node_state_prob(node_state_prob_init) {}


  bool is_invalid_prob;
  prob_t* tmp_store;
  prob_t** const branch_down_prob_storage;
  const time_gtor& tgm;
  const unsigned n_states;
  const prob_t* const * const leaf_distro;
  const prob_t* const * const branch_tprob;
  simple_array<prob_t> tprob_tmp;
  prob_t** node_state_prob;
};



// non-recursive portion of branch_parent distro evaluation:
//
static
void
get_parent_down_prob(const prob_t* transp_matrix,
                     const prob_t* child_distro,
                     prob_t* parent_distro,
                     tbdp_info& ti){

  linear_solver(transp_matrix,child_distro,parent_distro,ti.n_states);

  bool is_branch_invalid_prob(false);

  {
    bool is_nonzero_prob(false);
    for(unsigned i(0);i<ti.n_states;++i){
      if(parent_distro[i]<0.) {
        //        log_os << "NPROB: label,i,p(i): " << b->label << " " << i << " " << parent_distro[i] << "\n";

        // the current prob transition matrix could not have generated
        // the distro observed at the base node of this branch:
        parent_distro[i] = 0.;
        is_branch_invalid_prob=true;
      } else if(parent_distro[i]>0.) {
        is_nonzero_prob=true;
      }
    }

    if(! is_nonzero_prob) is_branch_invalid_prob=true;

    {
      prob_t sum(0.);
      for(unsigned i(0);i<ti.n_states;++i){ sum += parent_distro[i]; }
      if( std::fabs(sum -1.) > 1e-3 ) is_branch_invalid_prob = true;
    }

    if(is_nonzero_prob){
      pdistro_norm(parent_distro,parent_distro+ti.n_states);
    } else {
      pdistro_unif(parent_distro,parent_distro+ti.n_states);
    }
  }

  if(is_branch_invalid_prob) {
#ifdef USE_NNLS_ROOT
    for(unsigned ii(0);ii<100;++ii){
    int mda(static_cast<int>(ti.n_states));
    int m(mda);
    int n(mda);
    simple_array<double> b(m);
    simple_array<double> w(n);
    simple_array<double> zz(m);
    simple_array<double> A(m*m);
    simple_array<int> indx(n);
    simple_array<double> A(n*n);
    double err;
    int mode;

    std::copy(child_distro,child_distro+m,b.val);
    std::copy(transp_matrix,transp_matrix+m*m,A.val);

    nnls(A.val,mda,m,n,b.val,parent_distro,&err,w.val,zz.val,indx.val,&mode);

    if(mode != 1) { die("nnls failed"); }

    double sum(0.);
    for(unsigned i(0);i<n;++i){ sum += parent_distro[i]; }
      log_os << "sum: " << sum << "\n";
      log_os << "nnls sum,err,cerr " << sum << " " << err << " " << (err*err)/n << "\n";
    }

#else
    // linear solver failed to produce a valid probability at the
    // parent node, so iteratively find the parent distribution which
    // minimizes the child node least squared probability error
    //
    matrix_transpose_copy(transp_matrix,ti.tprob_tmp.ptr(),ti.n_states);
    //          pdistro_unif(parent_distro,ti.n_states);
#ifndef USE_CONJ_DIR_ROOT
    {
      for (unsigned ii(0);ii<100;++ii){
        parent_prob_leastsq_minfunc_indy ppm(ti.tprob_tmp.ptr(),child_distro,ti.n_states);

        const unsigned ndim(ppm.dim());

        simple_array<smlfloat> p(ndim);

        ppm.parent_prob_to_p(parent_distro,p.ptr());

        {
          static const smlfloat ntol(1e-28);
          smlfloat fmin;
          unsigned iter;
          smlfloat final_iter_delta_f;
          static const unsigned max_iter(5000);
          minimize_conj_gradient(p.ptr(),ppm,ntol,fmin,iter,final_iter_delta_f,max_iter);
#if 0
          double sum(0.);
          for(unsigned i(0);i<ndim;++i){ sum += parent_distro[i]; }
          log_os << "sum: " << sum << "\n"
                 << std::setprecision(20)
                 << "ROOT CONJ GRAD MINIMIZER ITERATIONS,fmin: " << iter << " " << fmin << "\n";
#endif
        }

        ppm.p_to_parent_prob(p.ptr(),parent_distro);
      }
    }
#else
    {
      parent_prob_leastsq_minfunc ppm(ti.tprob_tmp.val,child_distro,ti.n_states);

      const unsigned ndim(ppm.dim());

      simple_array<smlfloat> p(ndim);

      ppm.parent_prob_to_p(parent_distro,p.val);

      {
        static const smlfloat line_tol(1e-9);
        static const smlfloat start_ratio(0.05);
        static const smlfloat min_start_dist(1e-6);

        const smlfloat ntol(1e-28);
        smlfloat start_ntol(ntol);
        unsigned iter;
        smlfloat fmin;
        const unsigned max_iter(1000);

        smlfloat* conj_dir(new smlfloat[ndim*ndim]);
        std::fill(conj_dir,conj_dir+ndim*ndim,0.);
        for(unsigned i(0);i<ndim;++i) {
          const smlfloat start_dist( std::max(std::fabs(p.val[i]*start_ratio),min_start_dist) );
          conj_dir[i*(ndim+1)] = start_dist;
        }
        minimize_conj_direction(p.val,conj_dir,ppm,start_ntol,ntol,line_tol,fmin,iter,max_iter);
        log_os << std::setprecision(20);
        log_os << "CONJ DIR MINIMIZER ITERATIONS, fmin: " << std::min(iter+1,max_iter) << " " << fmin << "\n";
        delete [] conj_dir;
      }

      ppm.p_to_parent_prob(p.val,parent_distro);
    }
#endif
#endif
  }

  ti.is_invalid_prob = ti.is_invalid_prob || is_branch_invalid_prob;


#if DEBUG_BGR
  smlfloat sum(0.);
  for(unsigned i(0);i<ti.n_states;++i){
    sum += parent_distro[i];
  }
  log_os << "sum " << sum << "\n";

  log_os << "PUGGY " << b->label;
  for(unsigned i(0);i<3;++i){
    log_os << " " << parent_distro[i];
  }
  log_os << "\n";
#endif
}



// store the P(parent(state)|node,parent(node)) for each non-root node in tree
//
static
void
tree_branch_down_prob(const bi_tree_node* const b,
                      const prob_t*& branch_down_prob,
                      smlfloat& branch_weight,
                      tbdp_info& ti){

  const unsigned branch_id(b->branch_id());
  prob_t* branch_down_prob_tmp(ti.branch_down_prob_storage[branch_id]);
  branch_down_prob = branch_down_prob_tmp;

  const prob_t* child_distro(0);

  if(b==0){
    die("tree_branch_down_prob: null tree node");

  } else if(b->is_leaf()){  // leaf node
    branch_weight=1;
    const unsigned leaf_id(b->leaf_id());
    child_distro=ti.leaf_distro[leaf_id];

  } else { // interior node
    const prob_t* child1_prob(0);
    const prob_t* child2_prob(0);
    smlfloat child1_bw;
    smlfloat child2_bw;

    tree_branch_down_prob(b->child1(),child1_prob,child1_bw,ti);
    tree_branch_down_prob(b->child2(),child2_prob,child2_bw,ti);

    branch_weight=child1_bw+child2_bw;
    for(unsigned i(0);i<ti.n_states;++i) {
      ti.tmp_store[i] = (child1_prob[i]*child1_bw+child2_prob[i]*child2_bw)/branch_weight;
    }
    child_distro=ti.tmp_store;
  }

  if(ti.node_state_prob){
    for(unsigned i(0);i<ti.n_states;++i){
      ti.node_state_prob[branch_id+1][i] = child_distro[i];
    }
  }
  get_parent_down_prob(ti.branch_tprob[branch_id],child_distro,branch_down_prob_tmp,ti);
}



static
bool
subs_ml_model_root_node_prob_branch_tprob(const subs_ml_model& mdl,
                                          const prob_t* const * branch_tprob,
                                          prob_t* root_node_prob,
                                          const unsigned cat_no,
                                          const bool is_use_submodels,
                                          const unsigned submodel_no,
                                          prob_t** node_state_prob=0){

  const bi_tree& tree(mdl.tree());
  const time_gtor& tgm(mdl.get_time_gtor());
  const unsigned n_states(is_use_submodels ?
                          mdl.submodel_state_size(submodel_no) : mdl.state_size());
  const unsigned n_leaves(tree.leaf_size());
  const unsigned n_branches(tree.branch_size());
  const prob_t* const * leaf_distro(mdl.obs_seq_state_distro_cat(cat_no));

  simple_init_matrix<prob_t> tmp_leaf;
  if(is_use_submodels){
    tmp_leaf.init(n_leaves,n_states);
    for(unsigned i(0);i<n_leaves;++i){
      mdl.get_rate_gtor().submodel_pdistro_reduction(leaf_distro[i],submodel_no,tmp_leaf[i]);
    }
    leaf_distro=tmp_leaf.ptr();
  }

  if(tree.root()->is_leaf()) {
    const bool is_invalid_prob(false);
    for(unsigned i(0);i<n_states;++i) root_node_prob[i] = leaf_distro[0][i];
    return is_invalid_prob;
  }

  simple_matrix<prob_t> branch_down_prob_storage(n_branches,n_states);

  simple_array<prob_t> tmp_store(n_states);

  const bool init_is_invalid_prob(false);
  tbdp_info ti(init_is_invalid_prob,tmp_store.ptr(),branch_down_prob_storage.ptr(),
               tgm,n_states,leaf_distro,branch_tprob,node_state_prob);

  const prob_t* child1_prob(0);
  const prob_t* child2_prob(0);
  smlfloat child1_bw;
  smlfloat child2_bw;

  tree_branch_down_prob(tree.root()->child1(),child1_prob,child1_bw,ti);
  tree_branch_down_prob(tree.root()->child2(),child2_prob,child2_bw,ti);

  for(unsigned i(0);i<n_states;++i){ root_node_prob[i]=child1_prob[i]*child1_bw+child2_prob[i]*child2_bw; }
  pdistro_norm(root_node_prob,root_node_prob+n_states);

  if(ti.node_state_prob){
    for(unsigned i(0);i<ti.n_states;++i){
      ti.node_state_prob[0][i] = root_node_prob[i];
    }
  }

  if(ti.is_invalid_prob && mdl.opt().is_extra_warnings){
    warning("Iterative least-sq required to find root distribution");
  }

#ifdef DEBUG_BGR
  log_os << "PUPPY";
  for(unsigned i(0);i<3;++i){
    log_os << " " << root_node_prob[i];
  }
  log_os << "\n";
#endif

  return ti.is_invalid_prob;
}



static
bool
subs_ml_model_root_node_prob(const subs_ml_model& mdl,
                             prob_t* root_node_prob,
                             const unsigned cat_no,
                             const bool is_use_submodels,
                             const unsigned submodel_no
#ifdef BRANCH_SPECIFIC_OBS
                             ,
                             const prob_t * const * node_state_prob_in=0,
                             prob_t** node_state_prob_out=0
#endif
                             ){

#ifdef BRANCH_SPECIFIC_OBS
  workspace ws[3];
#else
  prob_t** node_state_prob_out(0);
#endif

  assert(cat_no<mdl.get_cat_manager().cat_size());

  const unsigned n_states(is_use_submodels ?
                          mdl.submodel_state_size(submodel_no) : mdl.state_size());
  const bi_tree& tree(mdl.tree());
  const unsigned n_branches(tree.branch_size());

  const unsigned N2(n_states*n_states);

  simple_matrix<prob_t> branch_tprob(n_branches,N2);

  const rates_func_options_base ropt(cat_no,false,false,is_use_submodels,submodel_no);

#ifdef BRANCH_SPECIFIC_OBS
  if(node_state_prob_in){
    site_model_edge_pdistros smp;
    smp.ropt_hookup(ropt);

    simple_array<prob_t> state_prob_tmp(n_states);
    simple_array<prob_t> tprob_tmp(n_states*n_states);
    simple_array<prob_t> tprob_tmp2(n_states*n_states);

    const SITE_MODEL::index_t sm(mdl.get_rate_gtor().site_model());
    for(unsigned b(0);b<n_branches;++b){
      static const unsigned n_sections(3);
      const smlfloat section_frac(1./static_cast<smlfloat>(n_sections));
      const smlfloat section_time(tgm.branch_time(b)*section_frac);

      for(unsigned j(0);j<n_sections;++j){
        const unsigned node_id(b+1);
        const unsigned parent_node_id(tree.get_parent_index(node_id));
        const smlfloat node_frac((section_frac*0.5)*(1+j*2));

        for(unsigned i(0);i<n_states;++i){
          state_prob_tmp.val[i] = node_state_prob_in[node_id][i]*(node_frac)+
            node_state_prob_in[parent_node_id][i]*(1.-node_frac);
        }
        set_site_model_edge_pdistros(sm,state_prob_tmp.val,smp);
        get_single_branch_subs_prob(tprob_tmp.val,section_time,mdl,ropt,ws);

        if(j==0){
          matrix_copy(tprob_tmp.val,branch_tprob.val[b],n_states);
        } else {
          matrix_mult(tprob_tmp.val,tprob_tmp2.val,branch_tprob.val[b],n_states);
        }
        if((j+1)!=n_sections){
          matrix_copy(branch_tprob.val[b],tprob_tmp2.val,n_states);
        }
      }
    }
  }
#else
  get_all_branch_subs_prob(branch_tprob.ptr(),mdl,ropt);
#endif


#ifdef DEBUG
  for(unsigned i(0);i<n_branches;++i){
    for(unsigned j(0);j<n_states;++j){
      pdistro_check(branch_tprob[i]+j*n_states,n_states,SUBS_ML_PTOL);
    }
  }
#endif

  const bool is_invalid_prob =
    subs_ml_model_root_node_prob_branch_tprob(mdl,branch_tprob.ptr(),root_node_prob,
                                              cat_no,is_use_submodels,submodel_no,
                                              node_state_prob_out);

#ifdef DEBUG
  pdistro_check(root_node_prob,n_states,SUBS_ML_PTOL);
#endif

  // shouldn't be required, but accumulated noise might be substantial with many cats, etc..
  pdistro_norm(root_node_prob,root_node_prob+n_states);

  return is_invalid_prob;
}


#ifdef BRANCH_SPECIFIC_OBS
bool
subs_ml_model_root_node_prob_nspround(const subs_ml_model& mdl,
                                      prob_t* root_node_prob,
                                      const unsigned cat_no,
                                      prob_t** node_state_prob){

  const unsigned n_submodels(mdl.submodel_size());
  const bool is_use_submodels(n_submodels>1);

  if(is_use_submodels){
   die("cannot collect node_state_probs in submodel mode");
  } else {
    const unsigned n_states(mdl.state_size());
    const unsigned n_nodes(mdl.tree().node_size());
    simple_matrix<prob_t> simple_nsp(n_nodes,n_states);

    bool is_invalid_round(false);
    static const unsigned n_rounds(2);
    for(unsigned i(0);i<n_rounds;++i){
      const prob_t* const * input_nsp(0);
      prob_t** output_nsp(0);
      if(i!=0){
        input_nsp=simple_nsp.val;
        for(unsigned n(0);n<n_nodes;++n){
          for(unsigned j(0);j<n_states;++j){
            simple_nsp.val[n][j] = node_state_prob[n][j];
          }
        }
      }
      if((i+1)!=n_rounds){
        output_nsp=node_state_prob;
      }
      bool is_invalid = subs_ml_model_root_node_prob(mdl,root_node_prob,
                                                     cat_no,is_use_submodels,0,
                                                     input_nsp,output_nsp);
      is_invalid_round = is_invalid_round || is_invalid;
    }
    return is_invalid_round;
  }
}
#endif


// this function just takes care of the submodel case
//
// outmost loop finds root distro for each submodel type separately
// and joins them.
//
// ..this can be overridden to find one big root distro, but will make
// the root solve the limiting step in the likelihood function for c5
// state
//
bool
subs_ml_model_root_node_prob(const subs_ml_model& mdl,
                             prob_t* root_node_prob,
                             const unsigned cat_no){

  const unsigned n_submodels(mdl.submodel_size());
  const bool is_use_submodels(n_submodels>1);

  if(is_use_submodels){
    bool is_invalid_prob(false);

    simple_array<prob_t*> submodel_root_node_prob(n_submodels);

    for(unsigned s(0);s<n_submodels;++s){
      const unsigned n_states(mdl.submodel_state_size(s));

      submodel_root_node_prob[s] = new prob_t[n_states];

      const bool iip =
        subs_ml_model_root_node_prob(mdl,submodel_root_node_prob[s],
                                     cat_no,is_use_submodels,s);

      is_invalid_prob = is_invalid_prob || iip;
    }

    mdl.get_rate_gtor().submodel_fuse_pdistros(submodel_root_node_prob.ptr(),root_node_prob);

    for(unsigned s(0);s<n_submodels;++s){
      delete [] submodel_root_node_prob[s];
    }

    return is_invalid_prob;
  } else {
    return subs_ml_model_root_node_prob(mdl,root_node_prob,cat_no,
                                        is_use_submodels,0);
  }
}

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

// $Id: lhood_model.h 1027 2007-11-26 22:16:31Z ctsa $

/// \file

#ifndef __SITE_PROB_MAKER_H
#define __SITE_PROB_MAKER_H

#include <memory>

struct condition_func;
struct lhood_model_prep_cat_site_prob;
struct site_data_fastlup;
struct subs_ml_model;


struct site_prob_maker {

  virtual ~site_prob_maker() {}

  virtual
  void
  get_cat_site_prob(const subs_ml_model& mdl,
                    const unsigned cat_no,
                    lhood_model_prep_cat_site_prob& nlpc) const = 0;
};



struct site_prob_maker_full : public site_prob_maker {

  explicit
  site_prob_maker_full(const site_data_fastlup& full_sdf)
    : site_prob_maker(), _full_sdf(full_sdf) {}

  virtual
  void
  get_cat_site_prob(const subs_ml_model& mdl,
                    const unsigned cat_no,
                    lhood_model_prep_cat_site_prob& nlpc) const;

private:
  const site_data_fastlup& _full_sdf;
};



struct site_prob_maker_root : public site_prob_maker {

  site_prob_maker_root(const site_data_fastlup& full_sdf,
                       const simple_matrix3d<smlfloat>& root_spp,
                       const std::auto_ptr<condition_func>* cf)
    : site_prob_maker(), _full_sdf(full_sdf), _root_spp(root_spp), _cf(cf) {}

  virtual
  void
  get_cat_site_prob(const subs_ml_model& mdl,
                    const unsigned cat_no,
                    lhood_model_prep_cat_site_prob& nlpc) const;

private:
  const site_data_fastlup& _full_sdf;
  const simple_matrix3d<smlfloat>& _root_spp;
  const std::auto_ptr<condition_func>* _cf;
};



struct site_prob_maker_rootcat : public site_prob_maker {

  site_prob_maker_rootcat(const site_data_fastlup& full_sdf,
                          const simple_matrix<smlfloat>& root_spp,
                          const simple_matrix<smlfloat>& othercat_site_prob,
                          const unsigned cat_no,
                          const condition_func& cf)
    : site_prob_maker(), _full_sdf(full_sdf), _root_spp(root_spp),
      _osp(othercat_site_prob), _cat_no(cat_no), _cf(cf) {}

  virtual
  void
  get_cat_site_prob(const subs_ml_model& mdl,
                    const unsigned cat_no,
                     lhood_model_prep_cat_site_prob& nlpc) const;

private:
  const site_data_fastlup& _full_sdf;
  const simple_matrix<smlfloat>& _root_spp;
  const simple_matrix<smlfloat>& _osp;
  const unsigned _cat_no;
  const condition_func& _cf;
};

#endif

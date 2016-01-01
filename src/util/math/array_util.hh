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

// $Id: array_util.hh 905 2007-10-10 22:47:06Z ctsa $

/// \file



/// \brief scale float array a to make E[a] == target_expect
///
/// works for unnormalized probability distribution
///
/// T must equal a floating point type
///
template <typename ForwardIterator,
          typename InputIterator,
          typename T>
T
array_scale_expectation(ForwardIterator a,
                        const ForwardIterator a_end,
                        InputIterator a_prob,
                        const T target_expect){

  T expect(0.),total_prob(0.);
  for(ForwardIterator a2(a);a2!=a_end;++a2){
    expect += (*a2) * (*a_prob);
    total_prob += *a_prob;
    ++a_prob;
  }
  expect /= total_prob;
  const T scale(target_expect/expect);

  for(;a!=a_end;++a) { *a *= scale; }
  return scale;
}



/// \brief scale selected members of float array a to make E[a] ==
/// target_expect
///
/// works for unnormalized probability distribution
///
/// T must equal a floating point type
///
template <typename ForwardIterator1,
          typename InputIterator,
          typename ForwardIterator2,
          typename T>
T
array_scale_expectation_selected_only(ForwardIterator1 a,
                                      const ForwardIterator1 a_end,
                                      InputIterator a_prob,
                                      ForwardIterator2 a_selected,
                                      const T target_expect){
  {  // if none selected, adjust all...

    bool is_none_selected(true);
    ForwardIterator2 asel(a_selected);
    for(ForwardIterator1 a2(a);a2!=a_end;++a2){
      if(*asel) { is_none_selected=false; break; }
      ++asel;
    }
    if(is_none_selected) {
      return array_scale_expectation(a,a_end,a_prob,target_expect);
    }
  }

  T expect(0.),selected_expect(0.),total_prob(0.);

  ForwardIterator2 asel(a_selected);
  for(ForwardIterator1 a2(a);a2!=a_end;++a2){
    total_prob += *a_prob;
    const T val((*a2) * (*a_prob));
    expect += val;
    if(*asel) selected_expect += val;
    ++a_prob;
    ++asel;
  }
  expect /= total_prob;
  selected_expect /= total_prob;

  const T selected_target_expect(target_expect-expect+selected_expect);
  const T scale(selected_target_expect/selected_expect);

  for(;a!=a_end;++a) {
    if(*a_selected) { *a *= scale; }
    ++a_selected;
  }
  return scale;
}

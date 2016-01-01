// -*- mode: c++; indent-tabs-mode: nil; -*-
//

/// \file permute_array.h
/// \brief
///


#ifndef __PERMUTE_ARRAY_H
#define __PERMUTE_ARRAY_H


/// \brief copy an array of T to all permutations of array
///
/// in_array size=N
/// out_array size=N*(N!)
///
/// This method stores all permutations in an array, therefore it is only
/// appropriate for small N (<=8 or so). A non-recursive permutation
/// iterator would be preferable for the general case.
///
template <typename T>
void
permute_array(T* in_array,
              const unsigned N,
              T* out_array);


#include "permute_array.hpp"

#endif

// -*- mode: c++; indent-tabs-mode: nil; -*-
//

/// \file permute_array.hpp
/// \brief
///


#include <algorithm>


template <typename T>
void
permute_array(T* in_array,
              const unsigned N,
              T* out_array){
  permute_array_impl(in_array,N,out_array,0);
}


template <typename T>
void
permute_array_impl(T* in_array,
                   const unsigned N,
                   T*& out_array,
                   const unsigned start){
  if(start==(N-1)){
    for(unsigned i(0);i<N;++i){
      out_array[i] = in_array[i];
    }
    out_array += N;
    return;
  }
  for(unsigned i(start);i<N;++i){
    std::swap(in_array[start],in_array[i]);
    permute_array_impl(in_array,N,out_array,start+1);
    std::swap(in_array[start],in_array[i]);
  }
}


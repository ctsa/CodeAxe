
#include "util/math/matrix_util.h"

#include <iostream>

int
main(int argc,char* argv[]){

  double A[4] = { 1.1,.4,2.,.8 };
  double v[2] = { .3,.5 };
  double z[2];

  matrix_vector_mult_sum(A,v,z,2,2);

  std::cout << z[0] << " " << z[1] << "\n";
  
};

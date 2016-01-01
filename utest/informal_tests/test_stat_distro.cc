// g++ -Idefault_include/ -L$CLAPACK_LIB_PATH -DUSE_LAPACK -DDEBUG -ggdb test.cc stationary_pdistro.cc -llapack -lcblaswr -lcblas -latlas -lF77
//

#include "stationary_pdistro.h"


#include <iostream>
main(){
  smlfloat mat[] = { 
    		      -0.2, 0.1,0.1, 
                     0.01, -0.02, 0.01, 
		      0.4, 0.5,-0.9 };
  smlfloat foo[] = { 1.,0.,0.};

  get_stationary_pdistro_from_rates(foo,mat,3);

  for(unsigned i(0);i<3;++i){
    std::cout << foo[i] << "\n";
  }
}

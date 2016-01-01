

#include "bioseq_util_pdf_conversions.h"

#include <iostream>


main(){
  double f[NSC4::SIZE];
  for(unsigned i(0);i<NSC4::SIZE;++i){
    if(! (std::cin >> f[i])){
      std::cerr << "distro incomplete!";
      exit(0);
    }
  }

  double x[NSCODON::SIZE];
  nsc4_pdf_2_nscodon_pdf(f,x);
  for(unsigned i(0);i<NSCODON::SIZE;++i) std::cout << x[i] << "\n";

  double x2[DINUC::SIZE];
  nsc4_pdf_2_codon_boundary_dinuc_pdf(f,true,x2);
  for(unsigned i(0);i<DINUC::SIZE;++i) std::cout << x2[i] << "\n";
}

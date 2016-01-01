

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

  double x2[DINUC::SIZE];
  nsc4_pdf_2_codon_boundary_dinuc_pdf(f,true,x2);

  double f2[NSC4::SIZE];
  nscodon_pdf_dinuc_pdf_2_nsc4_pre_pdf(x,x2,f2);

  double x3[NSCODON::SIZE];
  nsc4_pdf_2_nscodon_pdf(f2,x3);

  double x4[DINUC::SIZE];
  nsc4_pdf_2_codon_boundary_dinuc_pdf(f2,true,x4);

  for(unsigned i(0);i<NSCODON::SIZE;++i){
    if(std::fabs(x[i]-x3[i])>0.0001){
      std::cerr << "pooch! nscodon i,x,x3: " << i << " " << x[i] << " " << x3[i] << "\n";
    }
  }

  for(unsigned i(0);i<DINUC::SIZE;++i){
    if(std::fabs(x2[i]-x4[i])>0.0001){
      std::cerr << "pooch! dinuc i,x2,x4: " << i << " " << x2[i] << " " << x4[i] << "\n";
    }
  }

  for(unsigned i(0);i<NSC4::SIZE;++i){
    std::cout << "f/f2,i: " << i << " " << f[i] << " " << f2[i] << "\n";
  }
}

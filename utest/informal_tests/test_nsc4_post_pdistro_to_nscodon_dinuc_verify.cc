

#include "../bioseq_util_pdf_conversions.h"

#include <iostream>


main(){
  double x[NSCODON::SIZE];
  double x2[DINUC::SIZE];
  for(unsigned i(0);i<NSCODON::SIZE;++i){
    if(! (std::cin >> x[i])){
      std::cerr << "distro incomplete!";
      exit(0);
    }
  }

  for(unsigned i(0);i<DINUC::SIZE;++i){
    if(! (std::cin >> x2[i])){
      std::cerr << "distro incomplete!";
      exit(0);
    }
  }


  double f2[NSC4::SIZE];
  nscodon_pdf_dinuc_pdf_2_nsc4_post_pdf(x,x2,f2);

  double x3[NSCODON::SIZE];
  nsc4_pdf_2_nscodon_pdf(f2,x3);

  double x4[DINUC::SIZE];
  nsc4_pdf_2_codon_boundary_dinuc_pdf(f2,false,x4);

  double f3[NSC4::SIZE];
  nscodon_pdf_dinuc_pdf_2_nsc4_post_pdf(x3,x4,f3);

  double x5[NSCODON::SIZE];
  nsc4_pdf_2_nscodon_pdf(f3,x5);

  double x6[DINUC::SIZE];
  nsc4_pdf_2_codon_boundary_dinuc_pdf(f3,false,x6);


  for(unsigned i(0);i<NSCODON::SIZE;++i){
    if(std::fabs(x5[i]-x3[i])>0.0001){
      std::cerr << "pooch! nscodon i,x5,x3: " << i << " " << x5[i] << " " << x3[i] << "\n";
    }
  }

  for(unsigned i(0);i<DINUC::SIZE;++i){
    if(std::fabs(x6[i]-x4[i])>0.0001){
      std::cerr << "pooch! dinuc i,x6,x4: " << i << " " << x6[i] << " " << x4[i] << "\n";
    }
  }
#if 0
  for(unsigned i(0);i<NSC4::SIZE;++i){
    std::cout << "f/f2,i: " << i << " " << f[i] << " " << f2[i] << "\n";
  }
#endif
}

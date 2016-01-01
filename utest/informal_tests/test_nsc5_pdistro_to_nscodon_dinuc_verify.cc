

#include "../bioseq_util_pdf_conversions.h"
#include "../subs_ml_types.h"

#include <iostream>


// manipulate joint prob to have assigned marginal probabilities
//
static
void
shape_joint_pdistro_margins(const prob_t* bg_pdf_pos1,
                            const prob_t* bg_pdf_pos2,
                            prob_t p_joint[NUC::SIZE][NUC::SIZE]){

  static int MAXITER(100);
  for(int iter(0);true;iter++){
    if(iter>MAXITER){
      std::cerr << "Error:: can't shape joint nuc distro\n";
      abort();
    }

    // p(1,2) to p(1|2):
    for(int i(0);i<NUC::SIZE;++i){
      prob_t sum(0.);
      for(int j(0);j<NUC::SIZE;++j){ sum += p_joint[j][i]; }
      for(int j(0);j<NUC::SIZE;++j){ p_joint[j][i] /= sum; }
    }

    // test p(1)
    bool p_1_converge(true);
    for(int j(0);j<NUC::SIZE;++j){
      prob_t sum(0.);
      for(int i(0);i<NUC::SIZE;++i){
        sum += p_joint[j][i]*bg_pdf_pos2[i];
      }
      if(std::abs(sum-bg_pdf_pos1[j])>0.001){
        p_1_converge=false;
        break;
      }
    }

    // p(1|2) to p(1,2)
    for(int i(0);i<NUC::SIZE;++i){
      for(int j(0);j<NUC::SIZE;++j){
        p_joint[i][j] *= bg_pdf_pos2[j];
      }
    }

    // p(1,2) to p(2|1)
    for(int i(0);i<NUC::SIZE;++i){
      prob_t sum(0.);
      for(int j(0);j<NUC::SIZE;++j){ sum += p_joint[i][j]; }
      for(int j(0);j<NUC::SIZE;++j){ p_joint[i][j] /= sum; }
    }

    // test p(2)
    bool p_2_converge(true);
    for(int i(0);i<NUC::SIZE;++i){
      prob_t sum(0.);
      for(int j(0);j<NUC::SIZE;++j){
        sum += p_joint[j][i]*bg_pdf_pos1[j];
      }
      if(std::abs(sum-bg_pdf_pos2[i])>0.001){
        p_2_converge=false;
        break;
      }
    }

    // p(2|1) to p(1,2)
    for(int i(0);i<NUC::SIZE;++i){
      for(int j(0);j<NUC::SIZE;++j){
        p_joint[i][j] *= bg_pdf_pos1[i];
      }
    }

    if(p_1_converge && p_2_converge) return;
  }
}




main(){
  double f[NSC4::SIZE];
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

  double x2a[NUC::SIZE][NUC::SIZE];
  for(unsigned i(0);i<DINUC::SIZE;++i){
    NUC::index_t nx,n0;
    nx = DINUC::decode_nx(i);
    n0 = DINUC::decode_n0(i);
    x2a[nx][n0] = x2[i];
  }
  
  double n0[NUC::SIZE];
  nscodon_pdf_2_nuc_pdf_pos(x,0,n0);
  double n2[NUC::SIZE];
  nscodon_pdf_2_nuc_pdf_pos(x,2,n2);

  shape_joint_pdistro_margins(n2,n0,x2a);

  for(unsigned i(0);i<DINUC::SIZE;++i){
    NUC::index_t nx,n0;
    nx = DINUC::decode_nx(i);
    n0 = DINUC::decode_n0(i);
    x2[i] = x2a[nx][n0];
  }


  double f2[NSC5::SIZE];
  nscodon_pdf_dinuc_pdf_2_nsc5_pdf(x,x2,f2);

  double x3[NSCODON::SIZE];
  nsc5_pdf_2_nscodon_pdf(f2,x3);

  double x4[DINUC::SIZE];
  nsc5_pdf_2_codon_boundary_dinuc_pdf(f2,x4);

  double f3[NSC5::SIZE];
  nscodon_pdf_dinuc_pdf_2_nsc5_pdf(x3,x4,f3);

  double x5[NSCODON::SIZE];
  nsc5_pdf_2_nscodon_pdf(f3,x5);

  double x6[DINUC::SIZE];
  nsc5_pdf_2_codon_boundary_dinuc_pdf(f3,x6);


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

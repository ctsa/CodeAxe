// -*- mode: c++; indent-tabs-mode: nil; -*-
//

#include "bioseq_util.h"
#include "prob_util.h"
#include "stationary_pdistro.h"

#include <vector>

// group_seq specifies continuous sequence ranges: sequence
// dependencies are reset at the boundary between each range
//
void
test_nsc4_c4pre_distro(const double* nsc4_distro){

  static const unsigned size(100000000);

  vector<unsigned> seq;
  seq.resize(size);

  // cdf used for idependent 4mer draws:
  prob_t nsc4_cdf[NSC4::SIZE];
  pdistro_to_cdf(root_distro,nsc4_cdf,NSC4::SIZE);

  // cdf used for draws dependent on the 4mer's extra nuc:
  prob_t* nscodon_cdf_by_nx[NUC::SIZE];
  {
    prob_t nsc4_pdf[NSC4::SIZE];
    for(unsigned i(0);i<NSC4::SIZE;++i){ nsc4_pdf[i] = root_distro[i]; }

    for(unsigned i(0);i<NUC::SIZE;++i){
      nscodon_cdf_by_nx[i] = nsc4_pdf+(i*NSCODON::SIZE);
      pdistro_norm(nscodon_cdf_by_nx[i],NSCODON::SIZE);
      pdistro_to_cdf_inplace(nscodon_cdf_by_nx[i],NSCODON::SIZE);
    }
  }


  {
    // simulation values for C4PRE:
    int seq_start(0);
    int seq_end(size);
    int seq_increment(1);
    unsigned dependent_nuc(2);

    NUC::index_t n[CODON::BASE_SIZE];
    for(int i(seq_start); i != seq_end; i += seq_increment){
      if(i==seq_start){
        // for C4PRE(C4POST): first (last) position in each group is
        // not dependent on nuc from the preceding (following) codon:
        //
        seq[i] = NSC4::decode_nscodon(random_cdf_variate(nsc4_cdf,NSC4::SIZE));
      } else {
        // for C4PRE(C4POST): all other positions are found as a
        // function of the dependent_nuc in the preceding (following)
        // codon:
        //
        NSCODON::decode(n,static_cast<NSCODON::index_t>(seq[i-seq_increment]));
        seq[i] = random_cdf_variate(nscodon_cdf_by_nx[n[dependent_nuc]],NSCODON::SIZE);
      }
    }
  }



  prob_t nsc4_pdf[NSC4::SIZE];
  for(unsigned i(0);i<NSC4::SIZE;++i){ nsc4_pdf[i] = root_distro[i]; }

  prob_t stat_pdf[NSC4::SIZE];

  for(int i(0);i<NSC4::SIZE;++i){ stat_pdf[i] = 0.; }

  // convert new seq to pdf
  NUC::index_t nuc[CODON::BASE_SIZE];
  for(int i(1);i<size;++i){
    NSCODON::decode(nuc,static_cast<NSCODON::index_t>(seq[i-1]));
    const NSCODON::index_t c =  static_cast<NSCODON::index_t>(seq[i]);
    NSC4::index_t f = NSC4::encode(nuc[2],c);
    stat_pdf[f] += 1.;
  }

  for(int i(0);i<NSC4::SIZE;++i){ stat_pdf[i] /= static_cast<prob_t>(size-1); }

  pdistro_check(stat_pdf,NSC4::SIZE,1.e-5);


  for(int i(0);i<NSC4::SIZE;++i){
    NUC::index_t n;
    NSCODON::index_t c;
    NSC4::decode(n,c,i);
    std::cerr << "f: i,norm,stat: " << i << " " << NUC::syms[n] << "+" << NSCODON::print(c) << " "
              << nsc4_pdf[i] << " " << stat_pdf[i] << "\n";
  }

  // nx/n2
  prob_t nsc4_nx[NUC::SIZE];
  prob_t stat_nx[NUC::SIZE];
  prob_t nsc4_n2[NUC::SIZE];
  prob_t stat_n2[NUC::SIZE];
  for(int i(0);i<NUC::SIZE;++i){
    nsc4_nx[i] = 0.;
    stat_nx[i] = 0.;
    nsc4_n2[i] = 0.;
    stat_n2[i] = 0.;
  }

  for(int i(0);i<NSC4::SIZE;++i){
    NUC::index_t n;
    NSCODON::index_t c;
    NSC4::decode(n,c,i);
    NSCODON::decode(nuc,c);
    nsc4_nx[n] += nsc4_pdf[i];
    stat_nx[n] += stat_pdf[i];
    nsc4_n2[nuc[2]] += nsc4_pdf[i];
    stat_n2[nuc[2]] += stat_pdf[i];
  }

  std::cerr << "nx/n2:\n";
  for(int i(0);i<NUC::SIZE;++i){
    std::cerr << "c4 nx/n2 stat nx/n2: " << NUC::syms[i] << " "
              << nsc4_nx[i] << " " << nsc4_n2[i] << " "
              << stat_nx[i] << " " << stat_n2[i] << "\n";
  }

  std::cerr << "p(c|n)\n";
  for(int i(0);i<NSC4::SIZE;++i){
    NUC::index_t n;
    NSCODON::index_t c;
    NSC4::decode(n,c,i);
    std::cerr << "f: i,norm,stat: " << i << " " << NUC::syms[n] << "+" << NSCODON::print(c) << " "
              << nsc4_pdf[i]/nsc4_nx[n] << " " << stat_pdf[i]/stat_nx[n] << "\n";
  }


  //joint
  prob_t nsc4_joint[NUC::SIZE][NUC::SIZE];
  prob_t stat_joint[NUC::SIZE][NUC::SIZE];

  for(int i(0);i<NUC::SIZE;++i){
    for(int j(0);j<NUC::SIZE;++j){
      nsc4_joint[i][j] = 0.;
      stat_joint[i][j] = 0.;
    }
  }
  for(int i(0);i<NSC4::SIZE;++i){
    NUC::index_t n;
    NSCODON::index_t c;
    NSC4::decode(n,c,i);
    NSCODON::decode(nuc,c);
    nsc4_joint[n][nuc[2]] += nsc4_pdf[i];
    stat_joint[n][nuc[2]] += stat_pdf[i];
  }

  std::cerr << "p(nx,n2)\n";
  for(int i(0);i<NUC::SIZE;++i){
    for(int j(0);j<NUC::SIZE;++j){
    std::cerr << "c4 nx,n2 stat nx,n2 stat nx*n2: "
              << NUC::syms[i] << "," << NUC::syms[j] << " "
              << nsc4_joint[i][j] << " "
              << stat_joint[i][j] << " "
              << stat_nx[i]*stat_n2[j] << "\n";
    }
  }

  std::cerr << "p(n0,n1|n,n2)\n";
  for(int i(0);i<NSC4::SIZE;++i){
    NUC::index_t n;
    NSCODON::index_t c;
    NSC4::decode(n,c,i);
    std::cerr << "f: i,norm,stat: " << i << " " << NUC::syms[n] << "+" << NSCODON::print(c) << " "
              << nsc4_pdf[i]/nsc4_joint[n][nuc[2]] << " " << stat_pdf[i]/stat_joint[n][nuc[2]] << "\n";
  }

  //cond
  prob_t nsc4_cond[NUC::SIZE*NUC::SIZE];
  for(int i(0);i<NUC::SIZE;++i){
    prob_t sum(0.);
    for(int j(0);j<NUC::SIZE;++j){ sum += nsc4_joint[i][j];}
    for(int j(0);j<NUC::SIZE;++j){
      nsc4_cond[j+NUC::SIZE*i] = nsc4_joint[i][j]/sum;
    }
  }

  prob_t nsc4_stationary[NUC::SIZE];
  get_stationary_pdistro_from_pmatrix(nsc4_stationary,nsc4_cond,NUC::SIZE);

  std::cerr << "nxc4_stat/stat:\n";
  for(int i(0);i<NUC::SIZE;++i){
    std::cerr << "c4_stat real_stat: "
              << NUC::syms[i] << " "
              << nsc4_stationary[i] << " "
              << stat_nx[i] << "\n";
  }
}


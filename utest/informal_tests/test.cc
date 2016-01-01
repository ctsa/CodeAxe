
#include "util/bioseq_util_pdf_conversions.h"


main() {

  double nuc_all_pos[] = { 0.2,  0.3,  0.1,  0.4,
                           0.25, 0.35, 0.3,  0.1,
			   0.2,  0.4,  0.15, 0.25 };

  double nuc_all_pos2[NUC::SIZE*CODON::BASE_SIZE];
  double nscodon[NSCODON::SIZE];

  for(unsigned i(0);i<NUC::SIZE*CODON::BASE_SIZE;++i){
    nuc_all_pos2[i] = nuc_all_pos[i];
  }

  nuc_pdf_all_pos_2_nscodon_pdf(nuc_all_pos,nscodon);
  nscodon_pdf_2_nuc_pdf_all_pos(nscodon,nuc_all_pos);
  nuc_pdf_all_pos_2_nscodon_pdf(nuc_all_pos,nscodon);
  nscodon_pdf_2_nuc_pdf_all_pos(nscodon,nuc_all_pos2);

  for(unsigned i(0);i<NUC::SIZE*CODON::BASE_SIZE;++i){
    std::cerr << i << " " << nuc_all_pos[i] << " " << nuc_all_pos2[i] << "\n";
  }

}


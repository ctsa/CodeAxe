
#include "simple_util.h"
#include "util/math/prob_util.h"
#include "util/math/random_util.h"


#include <algorithm>
#include <iostream>


static
void
test_random_scaled_pdistro_variate_cached(){
  static const unsigned test_size(2000000);

  simple_array<double> test_pdistro(test_size);

  for(unsigned i(0);i<test_size;++i){
    double r;
    while((r=random_uniform()) <= 0.);
    test_pdistro.val[i] = r;
  }

  pdistro_norm(test_pdistro.val,test_pdistro.val+test_size);

  simple_array<double> draw_pdistro(test_size);
  

//static const unsigned mc[] = {64,4096,10000000};
  static const unsigned mc[] = {1024,10000000};
  static const unsigned n_mc(sizeof(mc)/sizeof(unsigned));

  for(unsigned j(0);j<n_mc;++j){
    random_scaled_pdistro_variate_cached<double> rspc(test_pdistro.val,test_size,mc[j]);

    std::fill(draw_pdistro.val,draw_pdistro.val+test_size,0.);
    for(unsigned i(0);i<test_size*2000000;++i){
      draw_pdistro.val[rspc.get_random_draw()] += 1.;
    }
    
    pdistro_norm(draw_pdistro.val,draw_pdistro.val+test_size);

    double sse(0.);
    double sae(0.);
    double max_ae(0.);
    unsigned max_ae_i(0);

    double lsse(0.);
    double lsae(0.);
    double max_lae(0.);

    for(unsigned i(0);i<test_size;++i){
      const double err(test_pdistro.val[i]-draw_pdistro.val[i]);
      sse += err*err;
      sae += std::abs(err);
      if(std::abs(err)>max_ae){
        max_ae=std::abs(err);
        max_ae_i=i;
      }

      static const double fudge(1e-12);
      const double lerr(std::log((draw_pdistro.val[i]+fudge)/(test_pdistro.val[i]+fudge)));
      lsse += lerr*lerr;
      lsae += std::abs(lerr);
      if(std::abs(lerr)>max_lae){
        max_lae=std::abs(lerr);
      }
    }
    sse /= static_cast<double>(test_size);
    sae /= static_cast<double>(test_size);
    lsse /= static_cast<double>(test_size);
    lsae /= static_cast<double>(test_size);

    std::cerr << "max_cdf,sse,rsse,sae,max_ae,max_ae_i: " << std::min(mc[j],test_size) << " " << sse << " " << std::sqrt(sse) 
              << " " << sae << " " << max_ae << " " << max_ae_i << "\n";
    std::cerr << "max_cdf,lsse,lrsse,lsae,max_lae: " << std::min(mc[j],test_size) << " " << lsse << " " << std::sqrt(lsse) 
              << " " << lsae << " " << max_lae << "\n";
  }
}


main(){
  long seed=random_init(1166125407);
  std::cerr << "seed: " << seed << "\n\n";
  test_random_scaled_pdistro_variate_cached();
}


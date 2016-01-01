
#include "chi2.h"

#include <iostream>

main(){
  std::cout << chi2(0.95,0.3) << "\n";
  std::cout << chi2(0.99,1) << "\n";
  std::cout << chi2(0.95,4) << "\n";
  std::cout << chi2(0.99,4) << "\n";
}



#include "cat_manager.h"


int main(){
  //expect success:
  cat_manager_optoins cmo;
  cat_manager x(cmo);

  //invalid cat_model: expect failure
  cmo.cat_model_str="0:gmm0,gmm1";
  cat_manager y(cmo);
}

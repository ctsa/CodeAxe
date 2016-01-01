
#include "bioseq_util.h"
#include "subs_ml_util.h"

#include <iostream>


main(){
  const unsigned site_id = get_site_id(NSC5::SIZE-1,NSC5::SIZE-1,NSC5::SIZE-1,NSC5::SIZE);
  std::cerr << "site_id " << site_id << "\n";

  unsigned s1,s2,o;
  decode_site_id(site_id,NSC5::SIZE,s1,s2,o);

  std::cerr << "s1,s2,o: " << s1 << " " << s2 << " " << o << "\n";
}

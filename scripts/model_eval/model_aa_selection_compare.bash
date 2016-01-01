#!/usr/bin/env bash

this_dir=$(dirname $0)
m2aa=$this_dir/model_to_aa_list.bash

m1=$1
m2=$2

join <($m2aa < $m1) <($m2aa < $m2) |\
awk '{ 
  diff=$2-$3; 
  ldiff=log($2)-log($3); 
  
  absdiff=diff; 
  if(absdiff<0) absdiff=-absdiff; 
    
  sum1+=absdiff/$2; 
  sum2+=(diff*diff); 
  sum3+=(ldiff*ldiff); 
  c+=1; 
} 

END { 
  print "absvar:",sum1/c,
        "var:",sum2/c,
        "stddev:",sqrt(sum2/c),
        "var(log(x)):",sum3/c,
        "stddev(log(x)):",sqrt(sum3/c), 
        "exp(stddev(log(x))):",exp(sqrt(sum3/c)); 
}'


#!/usr/bin/env bash

cat $HOME/proj/subs_ml/pos_select_test/test_aa_matrix/model.v2.mouse_rat_human.c4-pre.kimura.cpg-only.asymmetric.lockroot |\
grep aa_selection |\
awk '{$3="0"; print;}' >|\
tmp.aamat


awk '{
  if($2~/aa_selection/) {
    printf("rate_gtor_nscodon_base aa_selection_matrix_0_scale_param 1 1\n");
    while(getline < "tmp.aamat") print;
  }
  else {
    if($2=="select_model_0") $3=3;
    print;
  }
}'

rm tmp.aamat

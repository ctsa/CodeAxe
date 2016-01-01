#!/usr/bin/env bash



awk '{
  if($2~/aa_selection/){
    if(! aas_begin) { aas_begin=1; }
  } else {
    if(aas_begin){
      printf("rate_gtor_nscodon_base aa_selection_matrix_1_scale_param 0 1\n");
      printf("rate_gtor_nscodon_base aa_selection_matrix_2_scale_param 1 5.\n");
      aas_begin=0;
    }
  }

  if($2=="select_model_0") {
                         $3=3; print;
    $2="select_model_1"; $3=1; print;
    $2="select_model_2"; $3=1;
  }

  if($2=="n_site_select_matrix_cats") $3=3;
  if($2=="site_select_matrix_cat_label_1"){
    $3=0; print;
    $3=1; print;
    $3=2;
  }

  if($2=="site_select_matrix_cats_prior_1") {
    $3=1; $4=0.6; print;
    $3=1; $4=0.3; print;
    $3=1; $4=0.1;
  }

  if($2=="n_site_cats") $3=3;
  if($2=="site_cat_0_mr_ss_sm_ro") {
    $5=0; print;
    $5=1; print;
    $5=2;
  }

  print;
}'


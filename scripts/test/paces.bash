#!/usr/bin/env bash

# run through some basic routines
#

test_cmd(){
  cmd_str=$1
  if [ "$2" == "" ]; then 
    is_debug_bin=0
  else
    is_debug_bin=$2
  fi
  echo "testing: $cmd_str" >&2
  (  
    # if running a debug binary, then shut down stderr b/c there's just too much of it
    #
    if [ $is_debug_bin == 1 ]; then exec 2>&-; fi

    echo "$cmd_str" | bash
  )
}


mkdir -p test
file_tag=test/test


basic_test_track(){
  test_tag=$1
  model_flags="$2"
  data_cat_flags="$3"
  ml_flags="$4"
  bin=$5
  tree_file=$6
  sim_size=$7
  sim_type=$8

  tree_tag="-tree-file $tree_file"
  model_flags="$model_flags $tree_tag"

  model_tag=$file_tag.model.$test_tag
  data_tag=$file_tag.data.$test_tag

  data_sim=$data_tag.sim.gz
  out_data="-out >(gzip -c >| $data_sim)"
  in_data="-in-data <(gzip -dc $data_sim)"

  if [ "$sim_type" == "" ]; then sim_type=continuous; fi

  sim_options="-sim-assigned-cat -sim-model $sim_type"

  echo "test: $test_tag"
  if [ ! -e $model_tag.trained ]; then
    echo "starting..."
    test_cmd "$bin -create-model $model_flags -random-param >| $model_tag.random"
    test_cmd "$bin -norm -in-model $model_tag.random >| $model_tag.norm"
    test_cmd "$bin -sim $sim_options -sim-size $sim_size -in-model $model_tag.norm $out_data"
    test_cmd "$bin -create-model $model_flags $data_cat_flags $in_data >| $model_tag.start"
    test_cmd "$bin -create-model $model_flags $data_cat_flags $in_data -random-param >| $model_tag.start.r"
    
    is_debug_bin=0
    if [ $(nm $bin 2> /dev/null | wc -l) != 0 ]; then is_debug_bin=1; fi
 
    ml_group="-ml $ml_flags $in_data"
    test_cmd "$bin $ml_group -in-model $model_tag.start >| $model_tag.trained" $is_debug_bin
    test_cmd "$bin $ml_group -in-model $model_tag.start.r >| $model_tag.trained.r" $is_debug_bin
  fi

  for f in norm start trained start.r trained.r; do
    echo -n "${f}_lhood:	"
    test_cmd "$bin -lhood $ml_flags -in-model $model_tag.$f $in_data"
  done
}


bin_dir=../../bin
dbin=$bin_dir/CodeAxe.debug
xbin=$bin_dir/CodeAxe

tree_dir=$HOME/proj/subs_ml/all_kondra_runs/trees
tree_file1=$tree_dir/tree.human_chimp_mouse
tree_file2=$tree_dir/tree.primates

bigsize=20000
smallsize=5000

# simple nuc model
#
nuc_flags="-site-model nuc -rate-model k80" 
basic_test_track nucd "$nuc_flags" "" "-max-steps 20" $dbin $tree_file1 $smallsize iss

# nuc model with cats
#
cat_flags="-cat-model "0{smr{0}gmr{0}} 1{smr{1}gmr{0}} 2{smr{0}gmr{1}} 3{smr{1}gmr{1}}" -unlock-cat-prob"
basic_test_track nucd_cat "$nuc_flags $cat_flags" "" "-max-steps 20" $dbin $tree_file1 $smallsize discrete

# same nuc models with regular binary
#
basic_test_track nucx "$nuc_flags" "" "" $xbin $tree_file1 $bigsize iss
basic_test_track nucx_cat "$nuc_flags $cat_flags" "" "" $xbin $tree_file1 $bigsize discrete

# a nonrev nuc model w/ different roots:
#
nucnr_flags="-site-model nuc -rate-model nonrev"
basic_test_track nucxnr_root1 "$nucnr_flags -root-model obs-avg" "" "" $xbin $tree_file1 $bigsize iss
basic_test_track nucxnr_root2 "$nucnr_flags -root-model least-sq" "" "" $xbin $tree_file1 $bigsize iss
basic_test_track nucxnr_root3 "$nucnr_flags" "" "" $xbin $tree_file1 $bigsize iss
exit

# test new cat types:
#
tree_cat_flags1="-cat-model 0{sti{0}}1{sti{1}} -unlock-cat-prob"
tree_cat_flags2="-cat-model 0{gti{0}}1{gti{1}} -unlock-cat-prob"
global_cat_flags1="-cat-model 0{smo{0}}1{smo{1}} -unlock-cat-prob"
global_cat_flags2="-cat-model 0{gmo{0}}1{gmo{1}} -unlock-cat-prob"
basic_test_track nucx_tree_cat1 "$nuc_flags $tree_cat_flags1" "" "" $xbin $tree_file1 $bigsize discrete
basic_test_track nucx_tree_cat2 "$nuc_flags $tree_cat_flags2" "" "" $xbin $tree_file1 $bigsize discrete
basic_test_track nucx_global_cat1 "$nuc_flags $global_cat_flags1" "" "" $xbin $tree_file1 $bigsize discrete
basic_test_track nucx_global_cat2 "$nuc_flags $global_cat_flags2" "" "" $xbin $tree_file1 $bigsize discrete

# nuc model with cats - large tree
#
basic_test_track nucx_cat_bigtree "$nuc_flags $cat_flags" "" "" $xbin $tree_file2 $smallsize discrete

# simple dinuc model
#
dinuc_flags="-site-model dinuc -rate-model k80 -root-model obs-avg"
dinuc_ml_flags=""
basic_test_track dinuc "$dinuc_flags -context-model cpg-only" "" "$dinuc_ml_flags" $xbin $tree_file1 $bigsize discrete

# simple trinuc model
#
trinuc_flags="-site-model trinuc -rate-model k80 -root-model obs-avg"
basic_test_track trinuc "$trinuc_flags -context-model cpg-only" "" "$dinuc_ml_flags" $xbin $tree_file1 $bigsize discrete

# simple codon model
#
codon_flags="-site-model codon -rate-model k80 -select-model single -root-model obs-avg"
codon_ml_flags=""
basic_test_track codon "$codon_flags -context-model cpg-only" "" "$codon_ml_flags" $xbin $tree_file1 $bigsize discrete

# simple codon model w/ continuous simulator
#
basic_test_track codon_contsim "$codon_flags -context-model cpg-only" "" "$codon_ml_flags" $xbin $tree_file1 $bigsize

# simple codon model - large tree
#
basic_test_track codon_bigtree "$codon_flags -context-model cpg-only" "" "$codon_ml_flags" $xbin $tree_file2 $smallsize

# codon model with cats
#
codon_flags2="-site-model codon -rate-model gy94  -select-model symmetric -context-model indy -root-model obs-avg -reversible-tree"
basic_test_track codon_cat "$codon_flags2 $cat_flags" "" "$codon_ml_flags" $xbin $tree_file1 $bigsize discrete
basic_test_track codon_tree_cat1 "$codon_flags2 $tree_cat_flags1" "" "$codon_ml_flags" $xbin $tree_file1 $bigsize discrete
basic_test_track codon_tree_cat2 "$codon_flags2 $tree_cat_flags2" "" "$codon_ml_flags" $xbin $tree_file1 $bigsize discrete
basic_test_track codon_global_cat1 "$codon_flags2 $global_cat_flags1" "" "$codon_ml_flags" $xbin $tree_file1 $bigsize discrete

#codon_flags3="-site-model codon -rate-model k80  -select-model single -context-model indy"
basic_test_track codon_global_cat2 "$codon_flags2 $global_cat_flags2" "" "$codon_ml_flags" $xbin $tree_file1 $bigsize

exit

# codon model with assigned cats
#
basic_test_track codon_cat_acat "$codon_flags2 $tree_cat_flags" "-data-cat-map 0:1,1:2" "$codon_ml_flags" $xbin $tree_file1 $bigsize 0

# codon model with assigned cats and root categories
#
basic_test_track codon_cat_acat_roc "$codon_flags2 -site-select-matrix-cats 2 -group-select-cats 2 -site-select-matrix-cats-root-dup" "-site-select-matrix-cats-from-data -data-cat-map one:1,two:2" "$codon_ml_flags" $xbin $tree_file1 $bigsize 0

# c4 model with unassigned cats non-converging test 
#
c4_flags="-site-model c4-pre -rate-model nonrev  -select-model asymmetric -context-model triplet -root-model codon-dinuc"
c4_ml_flags="-max-steps 100"
basic_test_track c4pre "$c4_flags -site-select-cats 2 -group-rate-cats 2 -group-select-cats 2" "" "$c4_ml_flags" $xbin $tree_file1 $bigsize 1

# c5 model approx non-converging test 
#
c5_flags="-site-model c5 -rate-model nonrev  -select-model asymmetric -context-model factored-triplet -root-model codon-dinuc"
c5_ml_flags="-max-steps 50"
basic_test_track c5_approx "$c5_flags -c5-approx-model c4-gmean" "" "$c5_ml_flags" $xbin $tree_file1 $bigsize 1

# c5 model non-converging test 
#
c5_ml_flags="-max-steps 10"
basic_test_track c5 "$c5_flags" "" "$c5_ml_flags" $xbin $tree_file1 $bigsize 1

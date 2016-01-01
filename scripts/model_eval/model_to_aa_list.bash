#!/usr/bin/env bash
#
# extract model aa exchange parameters for [selection/substitution]
#
# usage: $0 [ -c ci_file ] [-l cat_label ] [ -s ]
#  -s = substitution rates (default: selection)
#

bin_dir=$(dirname $0)/../../bin
bin=$bin_dir/CodeAxe

is_subrate=0

while getopts ":c:l:s" Option
do
  case $Option in
  c ) ci_file=$OPTARG;;
  l ) cat_label=$OPTARG;;
  s ) is_subrate=1;;
  * ) echo "Unimplemented option chosen.";;   # DEFAULT
  esac
done



#optional ci file:
ciarg=""
is_ci=0
if [ "$ci_file" != "" ]; then
  ciarg="-in-confidence $ci_file"
  is_ci=1
fi


skey="SELECTION"
if [ $is_subrate == 1 ]; then skey="RATE"; fi

searchtag="SORTED AA $skey PARAMETERS"

$bin -report-model -in-model - $ciarg |\
awk -v cl=$cat_label -v stag="$searchtag" '
  BEGIN {
    if(cl != "") is_cl=1;
  }
  {
    if(scan) {
      if(!is_cl || cl_zone) {
        if(!NF) exit;
        print $0;
      } else {
        if(!NF) scan=0;
      }
    } 
    if(is_cl && $0~"SELECTION MATRIX CATEGORY:"){
      if($5 == cl) cl_zone=1;
      else         cl_zone=0;
    }
    if($0~stag) {
      getline;
      scan=1;
    }
  }' |\
tr '=' ' ' |\
sort |\
awk -v ic=$is_ci '{printf $1 " " $2; if(ic=="1") printf " " $4; printf ic; printf "\n";}'


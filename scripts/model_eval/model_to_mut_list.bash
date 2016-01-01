#!/usr/bin/env bash
#
# extract model cdm rates
#
# usage: $0 [ -c ci_file ] < model_file
#

bin_dir=$(dirname $0)/../../bin
bin=$bin_dir/CodeAxe

while getopts ":c:" Option
do
  case $Option in
  c ) ci_file=$OPTARG;;
  * ) echo "Unimplemented option chosen.";;   # DEFAULT
  esac
done



#optional ci file:
ciarg=""
if [ "$ci_file" != "" ]; then
  ciarg="-in-confidence $ci_file"
fi


skey="SELECTION"

searchtag="SORTED AA $skey PARAMETERS"

$bin -report-model -in-model - $ciarg |\
awk '/^CDM_RATE:/ {print $2,$3,$5}' |\
grep -v "N"


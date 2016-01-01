#!/usr/bin/env bash

my_dir=$(cd $(dirname $0); pwd)
bin=$my_dir/../../bin/CodeAxe

usage="Usage: $0 -o org_set -s site [ -a alignment -c category]"

align_type_base=clustal_rna_tuples

while getopts ":o:s:c:a:" Option
do
case $Option in
  o ) orgset_type=$OPTARG;;
  s ) site_type=$OPTARG;;
  c ) cat_type=$OPTARG;;
  a ) align_subtype=$OPTARG;;
  * ) echo $usage
      exit 1;;
esac
done


if [ $orgset_type"" == "" ] || [ $site_type"" == "" ]; then
  echo $usage
  exit 1
fi

if [ $align_subtype"" != "" ]; then
  align_subtype=.$align_subtype
fi

align_type=$align_type_base$align_subtype

orgset_tag=$(basename $orgset_type)
outfile=data.$orgset_tag$align_subtype.$site_type

datadir=$HOME/proj/subs_ml/genome_alignment/data
basedir=$datadir/$orgset_type
fsadirbase=$basedir/$align_type_base
fsadir=$basedir/$align_type

xtratags=""


if [ $cat_type"" != "" ]; then
  outfile=$outfile.$cat_type
  cattag=${fsadirbase}_cats/$cat_type
  xtratags="-in-cat $cattag -cat-labels $cattag/catlabel"
fi


cmd="$bin -site-model $site_type -process-seq -in-seq $fsadir $xtratags -out $outfile"
echo "****" $cmd
$cmd




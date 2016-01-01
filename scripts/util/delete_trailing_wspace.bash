#!/bin/bash

get_source() {
   ls {.,util,util/math}/*.{h,hh,hpp,cc,c} 2> /dev/null
   ls makefile 2> /dev/null 
}


for f in $(get_source); do
  cat $f |\
  sed 's/[ 	]*$//' >|\
  tmp

  if ! diff tmp $f > /dev/null; then 
    mv -f tmp $f
  else 
    rm tmp
  fi
done

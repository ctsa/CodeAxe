#!/usr/bin/env bash

get_source() {
  (cat $(echo {.,util/{bio,general,math}}/*.[hc]* |\
         tr ' ' '\n' |\
         awk '! /#/ {print}')) 2>| /dev/null
}

blank_filter() {
  awk '{if($0!="") print}'
}

simcomment_filter() {
  awk '{if($0!~/^ *\/\//) print}'
} 

if0_filter() {
awk '{ 
  if($0~/^ *#if/){
    l+=1;
    if(!bl){ if($0~/^ *#if 0/) {bl=l;} }
  }
  if(!bl) print $0;
  if($0~/^ *#else/){if(bl && (bl+1)>l) bl=0;}
  if($0~/^ *#endif/){l-=1; if(bl && bl>l) bl=0;}
}'
}

echo -n "naive sloc: "
get_source | wc -l


echo -n "nonblank/simplecomment/no_if0 sloc: "
get_source |\
blank_filter |\
simcomment_filter |\
if0_filter |\
wc -l


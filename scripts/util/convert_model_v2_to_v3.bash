#!/usr/bin/env bash

awk '\
{
  if($1=="subs_ml_model") { 
    if($2=="tree" || $2=="time") {
      $1="time_gtor"; 
      ++i; 
      buffer[i] = $0; 
    } else if($2=="rate_gtor_model") { 
      print $0; 
      for(j=1;j<=i;++j){ 
	print buffer[j];
      } 
    } else {
      print $0
    } 
  } else { 
    print $0; 
  } 
}'

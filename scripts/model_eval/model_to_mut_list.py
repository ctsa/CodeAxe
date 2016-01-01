#!/usr/bin/env python
"""
usage $0 model_file [ ci_file ] 
"""


import math
import sys

from model_compare_share import *



def get_ci_range(cifile) :

  ci_info = {}
  for line in open(cifile) :
    word = line.strip().split()
    if len(word) != 9 : continue
    if word[1] != "freeparam" : continue

    param_no=int(word[2])
    val=float(word[3])
    var=float(word[4])
    range=float(word[5])

    if var <= 0 : continue

    ci_info[param_no] = (val,range)

  return ci_info



def parse_mut_tag(tag) :
  w=tag.split('_')
  mut=w[5]
  fr=mut[0]
  to=mut[1:]
  ctext=w[7]
  if ctext == "indy" :
    label=fr+to
  else :
    label=ctext[0]+fr+ctext[2]+to
  return label



modelfile = sys.argv[1]

ci_info={}
if len(sys.argv) > 2 :
  cifile = sys.argv[2]
  ci_info=get_ci_range(cifile)
  

val_map=extract_model_free_param(modelfile)

k = val_map.keys()
k.sort()

for param_no in k :
  tag=val_map[param_no][0]
  val=val_map[param_no][1]

  ci_range=0
  if ci_info.has_key(param_no) :
    ci_range=ci_info[param_no][1]

  if tag.find("nuc_mutation_model") == -1 : continue

  label=parse_mut_tag(tag)

  sys.stdout.write("%s %f %f\n" % (label,val,ci_range))




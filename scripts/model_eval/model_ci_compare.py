#!/usr/bin/env python
"""
usage $0 reference_model_file trained_ci_file [ fake_trained_model_file ]
"""

import math
import sys

from model_compare_share import *


def get_ci_info(cifile) :

  ci_info = {}
  for line in open(cifile) :
    word = line.strip().split()
    if len(word) != 9 : continue
    if word[1] != "freeparam" : continue

    param_no=int(word[2])
    val=float(word[3])
    var=float(word[4])

    if var <= 0 : continue

    ci_info[param_no] = (val,var)

  return ci_info


if len(sys.argv) != 2 and len(sys.argv) != 3 :
  print __doc__
  sys.exit()

modelfile = sys.argv[1]
cifile = sys.argv[2]

is_fake_expect=False
if len(sys.argv)>3 :
  fake_expect_modelfile = sys.argv[3]
  is_fake_expect=True


ci_info=get_ci_info(cifile)

trueval_map=extract_model_free_param(modelfile)

fakeval_map={}
if is_fake_expect : fakeval_map=extract_model_free_param(fake_expect_modelfile)

k = trueval_map.keys()
k.sort()

for param_no in k :
  if ci_info.has_key(param_no) :
    tag=trueval_map[param_no][0]
    val=trueval_map[param_no][1]
    if is_fake_expect :
      train_val=fakeval_map[param_no][1]
    else :
      train_val=ci_info[param_no][0]
    train_var=ci_info[param_no][1]

    if train_var <= 0. : continue

    error=(train_val-val)/math.sqrt(train_var)
    sys.stdout.write("%i %s %f %f %f %f\n" % (param_no,tag,error,val,train_val,(val-train_val)))


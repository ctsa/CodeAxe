#!/usr/bin/env python
"""
usage $0 reference_model_file trained_model_file
"""

import math
import sys

from model_compare_share import *


modelfile = sys.argv[1]
modelfile2 = sys.argv[2]

trueval_map=extract_model_free_param(modelfile)
val_map=extract_model_free_param(modelfile2)

k = trueval_map.keys()
k.sort()

for param_no in k :
  if val_map.has_key(param_no) :
    tag=trueval_map[param_no][0]
    val=trueval_map[param_no][1]
    train_tag=val_map[param_no][0]
    train_val=val_map[param_no][1]

    if tag != train_tag : raise Exception("mismatched tags!")

    error="NA"
    sys.stdout.write("%i %s %s %f %f %f\n" % (param_no,tag,error,val,train_val,(val-train_val)))


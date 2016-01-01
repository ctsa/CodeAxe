#!/usr/bin/env python
"""
usage: $0 -model-file file -data-file file [-min-maxp arg1] < (alignment_id_no,XXX...) > (cat_no,alignment_id_no,XXX...)

  arg1 = value between 0. and 1. (default=0)
         this is the minimum posterior probability for a category to be called
"""

import os
import sys

from gpp_shared import *


def get_alignment_id(alignment_label) :
  fh=os.path.basename(alignment_label).split('.')[0]
  return int(fh[(fh.find("_")+1):])



def main() :
  model_file=data_file=""
  min_maxp_thresh=0.

  arg = sys.argv[1:]

  i = 0
  while i < len(arg)  :
    if   arg[i] == "-model-file" :
      i += 1
      model_file = arg[i]
    elif arg[i] == "-data-file" :
      i += 1
      data_file = arg[i]
    elif arg[i] == "-min-maxp" :
      i += 1
      min_maxp_thresh = float(arg[i])
    else :
      print __doc__
      return
    i += 1

  if model_file == "" or data_file == "" :
    print __doc__
    return
  if min_maxp_thresh < 0. or min_maxp_thresh > 1. :
    sys.exit("invalid min-maxp argument")

  fp=get_gpp_fp(model_file,data_file)

  catinfo=[]
  get_catinfo(fp,catinfo)

  id_to_cat={}

  for (label,max_cat) in get_gpp_max_cat_seq(fp,catinfo,min_maxp_thresh) :
    id = get_alignment_id(label)
    id_to_cat[id]=max_cat

  infp=sys.stdin
  outfp=sys.stdout

  for line in infp :
    id=int(line.strip().split()[0])
    if not id_to_cat.has_key(id) : 
      sys.stderr.write("id %i excluded from model\n" % id)
      continue
    outfp.write("%i %s" % (id_to_cat[id],line))


if __name__ == "__main__" : main()

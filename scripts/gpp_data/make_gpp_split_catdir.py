#!/usr/bin/env python
"""
usage: $0 -model-file file -data-file file -cat-dir dir [-min-maxp arg1]

  arg1 = value between 0. and 1. (default=0)
         this is the minimum posterior probability for a category to be called

new category info will be put into cat_dir
"""

import os
import sys

from gpp_shared import *


def main() :
  model_file=data_file=cat_dir=""
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
    elif arg[i] == "-cat-dir" :
      i += 1
      cat_dir = arg[i]
    elif arg[i] == "-min-maxp" :
      i += 1
      min_maxp_thresh = float(arg[i])
    else :
      print __doc__
      return
    i += 1

  if model_file == "" or data_file == "" or cat_dir == "" :
    print __doc__
    return

  if min_maxp_thresh < 0. or min_maxp_thresh > 1. :
    sys.exit("invalid min-maxp argument")

  fp=get_gpp_fp(model_file,data_file)

  catinfo=[]
  get_catinfo(fp,catinfo)

  if os.path.isdir(cat_dir) : raise Exception("cat_dir already exists: ",cat_dir)
  os.mkdir(cat_dir)

  for (label,max_cat) in get_gpp_max_cat_seq(fp,catinfo,min_maxp_thresh) :
    if max_cat == 0 : continue
  
    ofp=open(os.path.join(cat_dir,label+".cat"),"w")
    ofp.write("* %i\n" % max_cat)


  ncat=len(catinfo)

  clfp=open(os.path.join(cat_dir,"catlabel"),"w")
  clfp.write("# codon g"+str(ncat)+"mo group cat posterior prob parse\n#\n")
  for i,x in enumerate(catinfo) : clfp.write("%i %s\n" % (i+1,x[0]))


if __name__ == '__main__': main()


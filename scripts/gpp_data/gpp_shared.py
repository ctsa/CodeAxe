
import os
import sys


def get_catinfo(fp,catinfo) :
  "get category labels and prior probs"

  for line in fp :
    w=line.strip().split()
    if len(w) == 0 : continue

    if w[0] == "prior_prob:" :
      clabel=w[1]
      cpp=float(w[2])
      catinfo.append((clabel,cpp))
    if w[0] == "group_label" : break


def get_gpp_fp(model_file,data_file) :
  this_dir=os.path.dirname(sys.argv[0])
  cdir=os.path.join(this_dir,"../../bin")
  cbin=os.path.join(cdir,"CodeAxe")

  return os.popen(cbin+" -group-cat-post-prob -in-model "+model_file+" -in-data "+data_file,"r")


def get_gpp_max_cat_seq(fp,catinfo,min_maxp_thresh) :
  ncat=len(catinfo)
  cat_postp=[0.]*ncat

  for line in fp :
    w=line.strip().split()
    label=w[0]
    for i,p in enumerate(w[1:]) : cat_postp[i]=float(p)

    maxp=max(cat_postp)
    if maxp < min_maxp_thresh :
      max_cat=0
    else : 
      max_cat=cat_postp.index(maxp)+1

    yield (label,max_cat)


#!/usr/bin/env python
"""
usage $0 < mut_list > mut_list.decorated 
"""



def classify_mut_tag_ws(tag) :

  W=["A","T"]
  S=["C","G"]

  w=tag.split("->")

  to=w[1]
  if len(w[0]) == 1 :
    fr=w[0][0]
  else :
    fr=w[0][1]

  f="S"
  t="S"
  if fr in W : f="W"
  if to in W : t="W"

  return f+"->"+t



def classify_mut_tag_tiv(tag) :

  w=tag.split("->")

  to=w[1]
  if len(w[0]) == 1 :
    fr=w[0][0]
    context="X"+fr+"X"
  else:
    fr=w[0][1]
    context=w[0]

  is_cpg=(context.find("CG") != -1)

  is_ti=((fr == "A" and to == "G") or \
         (fr == "G" and to == "A") or \
         (fr == "C" and to == "T") or \
         (fr == "T" and to == "C"))

  if is_cpg :
    if is_ti : return "cpg_ti"
    else     : return "cpg_tv"
  else :
    if is_ti : return "ti"
    else     : return "tv"



import sys

infp=sys.stdin
outfp=sys.stdout

for line in infp :
  w=line.strip().split()
  if len(w) < 1 : continue

  wstag=classify_mut_tag_ws(w[0])
  titag=classify_mut_tag_tiv(w[0])
  
  outfp.write("%s %s %s" % (w[0],titag,wstag))
  if len(w) > 1 : outfp.write(" %s" % (" ".join(w[1:])))
  outfp.write("\n")


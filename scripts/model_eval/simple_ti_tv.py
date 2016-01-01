#!/usr/bin/env python
#


def add_rate(q,w) :
  if   (w[1].find("param_"+q[0]+"->"+q[1]+"_rate") != -1) or \
       (w[1].find("param_"+q[1]+"->"+q[0]+"_rate") != -1) : return float(w[3])
  else : return 0.


import sys

infp=sys.stdin
outfp=sys.stdout

ti=0.
tv=0.

for line in infp :
  if line.find("nuc_mutation_model") == -1 : continue
  w=line.strip().split()
  for q in [ ("A","G"),("C","T") ] :                     ti += add_rate(q,w)
  for q in [ ("A","C"),("A","T"),("C","G"),("G","T") ] : tv += add_rate(q,w)

outfp.write("transversion: %f\n" % (tv/8.))
outfp.write("transition: %f\n" % (ti/4.))


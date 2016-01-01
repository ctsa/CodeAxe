#!/usr/bin/env python
#
# $0 < model_file
# total branch time for each time-cat in model
#

import sys

bt={}

for line in sys.stdin :
  w=line.strip().split()
  if len(w) != 5 : continue

  if w[1][:9]=="time_cat_" :
    tgt=w[1][9:]
    tgt=tgt[:tgt.find("_branch")]
    if not bt.has_key(tgt) : bt[tgt] = 0.
    bt[tgt] += float(w[4])



kb=bt.keys()
kb.sort()

for k in kb : 
  sys.stdout.write("%s %f\n" % (k,bt[k]))


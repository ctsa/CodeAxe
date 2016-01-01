#!/usr/bin/env python

# convert 12*8 joint nuc params to 12*16 triplet nuc params:

import sys

infp = sys.stdin

paramno=0
p_joint=[0]*12*8
p_triplet=[0]*12*16

word = []

for line in infp :
  word = line.strip().split()
  if word[1].find("nuc_mutation_model") == -1 : continue
  param = float(word[3])

  p_joint[paramno] = param

  paramno += 1



for i in range (12*16) :
  base=i/16
  from0=i%4
  from2=(i/4)%4

  jointp1=base*8+from0
  jointp2=base*8+4+from2
  p_triplet[i] = p_joint[jointp1]*p_joint[jointp2]

  sys.stdout.write("%s xxx_%i 1 %f\n" % (word[0],i,p_triplet[i]))



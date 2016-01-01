#!/usr/bin/env python

import sys

infp=sys.stdin
outfp=sys.stdout

n = m = 0
sum1 = sum2 = sum3 = sum4 = 0.

for line in infp : 
  word=line.strip().split()

  true_val=float(word[3])
  estimate_val=float(word[4])

  abs_error=estimate_val-true_val
  if abs_error<0 : abs_error *= -1

  sum1 += abs_error
  sum2 += abs_error*abs_error
  n += 1

  if true_val > 0.0001 :
    abs_rel_error=abs_error/true_val
    rel_error2=abs_rel_error*abs_rel_error
    sum3+=abs_rel_error
    sum4+=rel_error2
    m += 1

import math

outfp.write("avg_abs_error: %f " % (sum1/float(n)))
outfp.write("error_variance: %f " % (sum2/float(n-1)))
outfp.write("error_stddev: %f " % math.sqrt(sum2/float(n-1)))
outfp.write("avg_abs_rel_error: %f " % (sum3/float(m)))
outfp.write("rel_error_variance: %f " % (sum4/float(m-1)))
outfp.write("rel_error_stddev: %f\n" % math.sqrt(sum4/float(m-1)))


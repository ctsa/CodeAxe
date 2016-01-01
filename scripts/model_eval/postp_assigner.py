#!/usr/bin/env python
"""
take site posteriour prob assignments and align them back to fasta sequence

usage: $0 postp fasta
"""


import sys

if len(sys.argv) != 3 :
  print __doc__
  sys.exit()

postp_file=sys.argv[1]
fasta_file=sys.argv[2]


pfp = open(postp_file)

site_size = int(pfp.readline().strip().split()[1])
site_repeat_size = int(pfp.readline().strip().split()[1])
site_start_offset = int(pfp.readline().strip().split()[1])

line = pfp.readline()
org_label = line.strip().split()

orglen = len(org_label)

site_lup = {}

for line in pfp :
  word=line.strip().split()
  sites=word[:orglen]
  cats=word[orglen:]

  # sort sites by org_label:
  tmp=zip(org_label,sites)
  tmp.sort()
  sites = [ x[1] for x in tmp ]
  sites = " ".join(sites)
  sites = sites.replace('+','')
  site_lup[sites] = cats



ffp = open(fasta_file)

seq=""
org_tag=""
seq_list = []

for line in ffp :
  if line[0] == '>' :
    if seq != "" : seq_list.append((org_tag,seq.upper()))
    header = line
    org_tag = line[1:].strip().split()[1]
    seq = ""
  else :
    seq += line.strip()

if seq != "" : seq_list.append((org_tag,seq.upper()))

seqlen=len(seq_list[0][1])



# first reorder seq_list to match sorted org_label order:
seq_list.sort()

miss=0

i=site_start_offset
n=0
while True :
  if (i+site_size)>seqlen : break

  if(i>=0) :
    sites = []
    for s in seq_list :
      sites.append(s[1][i:i+site_size])

    # if any N's are found in site, then it is interpreted as all 'N'
    for j in range(len(sites)) :
      is_ambig = False
      for s in sites[j] :
	if s == 'N' or s == '-' :
	  is_ambig = True

      if is_ambig :
	sites[j] = 'N'*len(sites[0])

    sites = " ".join(sites)
    if not site_lup.has_key(sites) :
      # ignore this at first and last site:
      if n != 0 and (i+site_repeat_size+site_size <= seqlen) :
        print "big trouble!",i,sites,seq_list[0]
        sys.exit()
    else :
      sys.stdout.write("%i %s\n" % (n+1," ".join(site_lup[sites])))

  i += site_repeat_size
  n += 1


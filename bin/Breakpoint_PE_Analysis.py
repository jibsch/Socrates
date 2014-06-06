#!/usr/bin/env python

import sys
import pysam
import math

if len(sys.argv) < 4:
  print "Usage: <bam file> <socrates output> <true mean>"
  sys.exit(0)

f=open(sys.argv[2])
soc_lines = f.readlines()
f.close()

MEAN = int(sys.argv[3])

def split_anchor(anchor):
  return [anchor.split(":")[0], int(anchor.split(":")[1])]
def mean(inserts):
  m = 0
  for i in inserts:
    m += i
  return m/len(inserts)
def std(insert, mean):
  s = 0
  for i in inserts:
    s += math.pow(i-mean, 2)
  return math.sqrt( s/len(inserts) )
def t_test(sample_mean, true_mean, sample_size, sample_std):
  return (sample_mean - true_mean) / (sample_std * math.sqrt(sample_size) )

samfile = pysam.Samfile(sys.argv[1], 'rb')

if soc_lines[0].startswith("#"):
  print soc_lines[0]
  soc_lines = soc_lines[1:]

for soc in soc_lines:
  anchor1 = split_anchor(soc.split('\t')[0])
  anchor1.append(soc.split('\t')[1])
  anchor2 = split_anchor(soc.split('\t')[3])
  anchor2.append(soc.split('\t')[4])

  anchor = anchor1
  #short_sc_cluster = False
  out = []
  for i in [1,2]:
    #print anchor 
    inserts = []
    if anchor[2] == '-':
      true_mean = MEAN
      for r in samfile.fetch(anchor[0], anchor[1]-1, anchor[1]) :
        s = r.pos +1
        if s == anchor[1] and  r.isize > 0 and not r.is_reverse:
          inserts.append(r.isize)

    else :
      true_mean = -MEAN
      for r in samfile.fetch(anchor[0], anchor[1]-1, anchor[1]) :
        e = r.aend
        if e == anchor[1] and  r.isize < 0 and r.is_reverse:
          #inserts.append(r.isize )
          isize = r.isize - (100 - r.qlen) + r.qstart
          inserts.append(isize)

    if len(inserts) != 0:
      sample_mean =  mean(inserts)
      sample_std = std( inserts, sample_mean)
      if sample_std == 0:
        t = 0
      else :
        t = t_test(sample_mean, true_mean, len(inserts), sample_std)
      out.append("mean=%f, std=%f, t=%f" %(sample_mean, sample_std, t) )
    else : out.append( "X" )

    anchor = anchor2
   # short_sc_cluster = (soc.split('\t')[24] == 'short SC cluster\n')
    if soc.split('\t')[24] == 'short SC cluster\n' :
      out.append( "X" )
      break

  print out

samfile.close()


#!/usr/bin/env python

import bisect
import sys

def close(pos, exon, bp_direction, margin):
  if pos >= exon[0] - margin and pos <= exon[0] + margin and bp_direction == "-":
    return True
  elif pos >= exon[1] - margin and pos <= exon[1] + margin and bp_direction == "+":
    return True
  else:
    return False

if len(sys.argv) < 2: 
  print "Usage: breakpoint coordinate file please!"
  sys.exit(0)


f=open('/home/users/allstaff/schroeder/references/human/refGene.txt')
lines=f.readlines()
f.close()

exons={}

for l in lines:
  sp=l.split()
  name=sp[1]
  chr=sp[2]
  if not chr in exons:
    exons[chr] = []
  ex = exons[chr]
  st=sp[9].split(",")
  en=sp[10].split(",")
  exon_no = 1
  for s,e in zip(st[:-1],en):
    bisect.insort_right( ex, (int(s), int(e), name, exon_no) )
    exon_no+=1

f=open(sys.argv[1])
lines=f.readlines()
f.close()
  
exex = 0
exon_boundary_margin = 6

for l in lines :
  sp=l.split()
  c1=sp[0]
  p1=int(sp[1])
  d1=sp[2]
  c2=sp[3]
  p2=int(sp[4])
  d2=sp[5]
  if  c1 not in exons or c2 not in exons:
    print "no gene connection for line: ",
    print l,
    continue
  i = bisect.bisect_left(exons[c1], (p1,p1,'',1))
  j = 1
  genes1 = []
  while i-j >= 0 and close(p1, exons[c1][i-j], d1, exon_boundary_margin ):
     genes1.append(exons[c1][i-j])
     j += 1
  j = 0
  while i+j < len (exons[c1]) and close(p1, exons[c1][i+j], d1, exon_boundary_margin ):
    genes1.append(exons[c1][i+j])
    j += 1



  i = bisect.bisect_left(exons[c2], (p2,p2,'',1))
  j = 1
  genes2 = []
  while i-j >= 0 and close(p2, exons[c2][i-j], d2,exon_boundary_margin  ):
     genes2.append(exons[c2][i-j])
     j += 1
  j = 0
  while i+j < len (exons[c2]) and close(p2, exons[c2][i+j], d2, exon_boundary_margin ):
    genes2.append(exons[c2][i+j])
    j += 1

  #print "genes for %d: %s\ngenes for %d: %s\n------" %(p1, str(genes1), p2, str(genes2))

  if len(genes1) == 0 and len(genes2) == 0:
    print "no gene connection for line: ",
    print l,
    continue
  if len(genes1) == 0:
    print "partial gene link into %s: " %(str(genes2)),
    print l,
    continue
  if len(genes2) == 0:
    print "partial gene link into %s: " %(str(genes1)),
    print l,
    continue

  #g2=[x[2] for x in genes2]
  g2 = {x[2]:x for x in genes2}
  link = []
  for x in genes1:
    if x[2] in g2:
      #print "breakpoint is exon exon link"
      y = g2[x[2]]
      link.append((x[2],x[3],y[3]))
  if link ==[]:
    print "breakpoint links two different genes: %s VS %s" %(str(genes1), str(genes2))
  else : 
    exex += 1
    print "exon exon junction: ",
    print link

#print ""
#print "Total exon junctions: %d" %(exex)

     

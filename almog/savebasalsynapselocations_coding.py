from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import time
import sys
import random

random.seed(1)

maxLen = 356
Nsyns = 1000

lenToStart = 0.0
lenToEnd = 356.0
maxSynsPerSeg = 17

synlocs = []

h("""
load_file("myrun.hoc")
objref cvode, tvec
cvode = new CVode()
cvode.active(1)
cvode.atol(0.001)

access a_soma
double siteVec[2]
tvec = new Vector()
objref sl,syns["""+str(Nsyns)+"""]
sl = new List()
""")
synsInSegs = [0]*len(h.dend)
for istim in range(0,Nsyns):
  myiseg = -1
  while myiseg == -1:
    x = lenToStart+(maxLen-lenToStart)*random.random()
    h("""sl = locateSites(\"dend\","""+str(x)+""")
Nsegs_x = sl.count()
""")
    iseg = random.randint(0,h.Nsegs_x-1)
    h("myseg = sl.o["+str(iseg)+"].x[0]")
    if synsInSegs[int(h.myseg)] < maxSynsPerSeg:
      myiseg = int(h.myseg)
      break
    print "istim = "+str(istim)+", x = "+str(x)+", continue searching for myseg..."
  synsInSegs[myiseg] = synsInSegs[myiseg] + 1
  h("""
siteVec[0] = sl.o["""+str(iseg)+"""].x[0]
siteVec[1] = sl.o["""+str(iseg)+"""].x[1]
dend[siteVec[0]] {
  syns["""+str(istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(istim)+"""].e = 0
  syns["""+str(istim)+"""].tau = 5
  syns["""+str(istim)+"""].onset = 10000
}
""")
  synlocs.append([h.siteVec[0],h.siteVec[1]])
  
picklelist = [Nsyns,maxSynsPerSeg,maxLen,synlocs]
file = open('basalsynlocs'+str(lenToStart)+'-'+str(lenToEnd)+'.sav', 'w')
pickle.dump(picklelist,file)
file.close()

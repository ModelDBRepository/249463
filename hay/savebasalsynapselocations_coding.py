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

proximalpoint = 400
distalpoint = 620
BACdt = 5.0
fs = 8

lenToStart = 0.0
lenToEnd = 282.0
maxSynsPerSeg = 50
Nsyns = 1000

maxSynsPerSeg = 12

synlocsAll = []

for icell in range(0,1):
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"
  synlocs = []

  h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
objref st1, st2
st1 = new IClamp(0.5)
st2 = new IClamp(0.5)
L5PC.soma st1
L5PC.soma st2
objref vsoma, vdend, recSite, vdend2, isoma, cadend, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
objref sl,ns,syn1,con1,isyn, tvec, syns["""+str(Nsyns)+"""]
isyn = new Vector()
tvec = new Vector()
sl = new List()
double siteVec[2]
}
""")
  synsInSegs = [0]*len(h.L5PC.dend)
  for istim in range(0,Nsyns):
    myiseg = -1
    while myiseg == -1:
      x = lenToStart+(lenToEnd-lenToStart)*random.random()
      h("""sl = L5PC.locateSites("dend","""+str(x)+""")
Nsegs_x = sl.count()
""")
      iseg = random.randint(0,h.Nsegs_x-1)
      h("iseg = sl.o["+str(iseg)+"].x[0]")
      if synsInSegs[int(h.iseg)] < maxSynsPerSeg:
        myiseg = int(h.iseg)
        break
      print "istim = "+str(istim)+", x = "+str(x)+", continue searching for iseg..."
    synsInSegs[myiseg] = synsInSegs[myiseg] + 1
    h("""
siteVec[0] = sl.o["""+str(iseg)+"""].x[0]
siteVec[1] = sl.o["""+str(iseg)+"""].x[1]
access L5PC.dend[siteVec[0]]
L5PC.dend[siteVec[0]] {
  syns["""+str(istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(istim)+"""].e = 0
  syns["""+str(istim)+"""].tau = 5
  syns["""+str(istim)+"""].onset = 10000 + """+str(BACdt)+""" 
}
""")
    synlocs.append([h.siteVec[0],h.siteVec[1]])
  
  h("""
access L5PC.soma
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
sl = new List()
access L5PC.soma
isoma = new Vector()
cvode.record(&st1.i,isoma,tvec)
""")

  synlocsAll.append(synlocs[:])
picklelist = [Nsyns,maxSynsPerSeg,[],synlocsAll]
file = open('basalsynlocs'+str(lenToStart)+'-'+str(lenToEnd)+'.sav', 'w')
pickle.dump(picklelist,file)
file.close()

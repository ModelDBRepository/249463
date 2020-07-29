import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import time
import sys
import random
from os.path import exists

random.seed(1)

v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
BACdt = 5.0
fs = 8
ITERS = 30
tstop = 11000.0

import mutation_stuff
MT = mutation_stuff.getMT()
defVals = mutation_stuff.getdefvals()
keyList = defVals.keys()
for idefval in range(0,len(keyList)):
  if type(defVals[keyList[idefval]]) is not list:
    defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments
unpicklefile = open('scalings_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
theseCoeffsAllAll = unpickledlist[0]
theseMutValsAllAll = unpickledlist[2]

lensToStart = [100.0 + x*50 for x in range(0,16)]

for istartdist in range(0,len(lensToStart)):
 startdist = lensToStart[istartdist]
 unpicklefile = open('synlocs'+str(startdist)+'.sav', 'r')
 unpickledlist = pickle.load(unpicklefile)
 unpicklefile.close()
 Nsyns = unpickledlist[0]
 synlocsAll = unpickledlist[3]
 startdist = int(startdist)

 maxLens = [1300,1185]

 gsAllAll = []

 for icell in range(0,7):
  synlocs = synlocsAll[0]
  gsAll = []

  theseCoeffsAll = theseCoeffsAllAll[icell]

  coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

  counter = -1
  
  for igene in range(0,len(MT)):
   gsThisGene = []
   for imut in range(0,len(MT[igene])):
    gsThisMut = []
    nVals = len(MT[igene][imut])*[0]
    thesemutvars = []
    theseCoeffs = theseCoeffsAll[igene][imut]
    for imutvar in range(0,len(MT[igene][imut])):
      thesemutvars.append(MT[igene][imut][imutvar][0])
      if type(MT[igene][imut][imutvar][1]) is int or type(MT[igene][imut][imutvar][1]) is float:
        MT[igene][imut][imutvar][1] = [MT[igene][imut][imutvar][1]]
      nVals[imutvar] = len(MT[igene][imut][imutvar][1])
    cumprodnVals = cumprod(nVals)
    allmutvars = cumprodnVals[len(MT[igene][imut])-1]*[thesemutvars]
    allmutvals = []
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      allmutvals.append([0]*len(thesemutvars))
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      for imutvar in range(0,len(MT[igene][imut])):
        if imutvar==0:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][iallmutval%nVals[imutvar]]
        else:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][(iallmutval/cumprodnVals[imutvar-1])%nVals[imutvar]]
      
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      counter = counter + 1                                                                                                                                                               

      gsThisMutVal = []

      if exists('thresholddistalamp'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav'):
        print 'thresholddistalamp'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav exists'
        unpicklefile = open('thresholddistalamp'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav', 'r')
        unpickledlist = pickle.load(unpicklefile)
        unpicklefile.close()
        gsThisMutVal = unpickledlist[1]
      #picklelist = [theseCoeffsAll,gsThisMutVal,MT]
      #file = open('thresholddistalamp'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav', 'w')
      #pickle.dump(picklelist,file)
      #file.close()
      gsThisMut.append(gsThisMutVal[:])
    gsThisGene.append(gsThisMut[:])
   gsAll.append(gsThisGene[:])
  gsAllAll.append(gsAll[:])
  
 picklelist = [theseCoeffsAllAll,gsAllAll,MT]
 file = open('thresholddistalamp'+str(startdist)+'.sav', 'w')
 pickle.dump(picklelist,file)
 file.close()

 #picklelist = [theseCoeffsAllAll,gsThisAllAll,MT]
 #file = open('thresholddistalamp'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav', 'w')
 #pickle.dump(picklelist,file)
 #file.close()
  

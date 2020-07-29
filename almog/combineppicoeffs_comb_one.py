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
from setparams import *
from os.path import exists

random.seed(1)
import mutation_stuff
MT = mutation_stuff.getMT()

unpicklefile = open('scalings_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()

theseCoeffsAllAll = unpickledlist[0]
theseMutValsAllAll = unpickledlist[2]
ISIs = [10.0*x for x in range(0,51)]
lensToStart = [150.0, 300.0, 450.0, 600.0]

icomb = 0
#combs_all = [ [[0,5,1], [1,2,15], [2,4,7], [3,1,0], [5,0,0], [8,5,0], [13,2,0]],           #max, Hay cell 0
#              [[0,5,0], [1,3,0], [2,5,1], [3,1,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],  #min, Hay cell 0
#              [[0,5,1], [1,2,15], [2,4,7], [3,1,0], [5,0,0], [8,5,0], [13,3,0]],           #max, Hay cell 1
#              [[0,5,0], [1,3,0], [2,5,1], [3,1,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],  #min, Hay cell 1
#              [[0,5,1], [1,2,15], [2,4,7], [3,1,0], [5,0,0], [8,5,0], [13,2,0]],           #max, Hay cell 2
#              [[0,5,0], [1,3,0], [2,5,1], [3,1,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],  #min, Hay cell 2
#              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,5,0]],           #max, Hay cell 3
#              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,4,0]],  #min, Hay cell 3
#              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,1,0]],           #max, Hay cell 4
#              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],  #min, Hay cell 4
#              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,5,0]],           #max, Hay cell 5
#              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,3,0]],  #min, Hay cell 5
#              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,0,0]],           #max, Hay cell 6
#              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]] ] #min, Hay cell 6                                                                                                  
combs_all = [ [[0,2,9], [1,2,9], [2,4,1], [3,1,0], [5,0,0], [8,5,0], [13,0,0]],           #max, Hay cell 0
              [[0,2,4], [1,2,8], [2,5,0], [3,0,0], [6,0,0], [8,2,0], [12,1,1], [13,5,0]], #min, Hay cell 0
              [[0,4,5], [1,2,3], [2,1,7], [3,0,1], [5,0,0], [8,2,0], [12,1,1]],           #max, Almog cell 0
              [[0,3,2], [1,2,14], [2,4,4], [3,1,0], [6,1,0], [8,5,0], [12,0,0], [13,5,0]] ]

if len(sys.argv) > 1:
  icomb = int(float(sys.argv[1]))
combs = combs_all[icomb]


for istartdist in range(0,len(lensToStart)):
 startdist = lensToStart[istartdist]
 startdist = int(startdist)
 if len(sys.argv) > 2 and int(float(sys.argv[2])) != istartdist:
   continue

 for icell in range(0,7):
  theseCoeffsAll = theseCoeffsAllAll[icell]
  gCoeffsAll = []
  if len(sys.argv) > 3 and int(float(sys.argv[3])) != icell:
    continue

  unpicklefile = open('thresholddistalamp'+str(startdist)+'_cs'+str(icell)+'_comb'+str(icomb)+'.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  gs_thiscomb = unpickledlist[1]

  iitercounter = -1
  gCoeffsThisMut = []
  for iter in [0, 2, 6, 8, -1]:
    iitercounter = iitercounter + 1
    gCoeffsThisIter = []
    thisg = gs_thiscomb[iitercounter][-1]                                                                                                                             
    counter = -1 

    if not exists('PPIcoeffs'+str(startdist)+'_cs'+str(icell)+'_comb'+str(icomb)+'_iter'+str(iter)+'.sav'):
      print 'PPIcoeffs'+str(startdist)+'_cs'+str(icell)+'_comb'+str(icomb)+'_iter'+str(iter)+'.sav does not exist, continuing'
      continue

    unpicklefile = open('PPIcoeffs'+str(startdist)+'_cs'+str(icell)+'_comb'+str(icomb)+'_iter'+str(iter)+'.sav', 'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    gCoeffsThisMutVal = unpickledlist[1][0]

    gCoeffsThisMut.append(gCoeffsThisMutVal[:])

  picklelist = [theseCoeffsAll,gCoeffsThisMutVal,ISIs,MT]
  file = open('PPIcoeffs'+str(startdist)+'_cs'+str(icell)+'_comb'+str(icomb)+'.sav', 'w')
  pickle.dump(picklelist,file)
  file.close()

#copied from calcifcurves2.py 30.10.2018
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import sys
from setparams import *
from os.path import exists
import time
from os.path import exists

v0 = -62
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
Is = [0.69+0.005*x for x in range(0,43)]
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

import mutation_stuff
MT = mutation_stuff.getMT()
defVals = mutation_stuff.getdefvals()
keyList = defVals.keys()
mySuffixes = mutation_stuff.getsuffixes()
mySuffixExceptions = mutation_stuff.getsuffixexceptions()

printcounter = 0;
print("printcounter="+str(printcounter)); printcounter = printcounter+1
unpicklefile = open('scalings_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()

theseCoeffsAllAll = unpickledlist[0]
theseMutValsAllAll = unpickledlist[2]

print("printcounter="+str(printcounter)); printcounter = printcounter+1
print("printcounter="+str(printcounter)); printcounter = printcounter+1

paramdicts = []
paramdicts.append({})                                                                                                                                         # 4-6 spikes per burst, control
paramdicts.append({'transvec.x(31)': 1.25, 'transvec.x(32)': 1.25})                                                                                           # 4-5 spikes per burst         
paramdicts.append({'transvec.x(31)': 1.5, 'transvec.x(32)': 1.5})                                                                                             # 3-4 spikes per burst         
paramdicts.append({'transvec.x(31)': 2.0, 'transvec.x(32)': 2.0})                                                                                             # 3-4 spikes per burst        
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0})                                                                                             # 2-3 spikes per burst       
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0, 'transvec.x(20)': 1.3, 'transvec.x(21)': 1.3, 'transvec.x(25)': 1.3, 'transvec.x(26)': 1.3}) # 2 spikes per burst        
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0, 'transvec.x(20)': 1.6, 'transvec.x(21)': 1.6, 'transvec.x(25)': 1.6, 'transvec.x(26)': 1.6}) # 1-2 spikes per burst     

icomb = 0
combs_all = [ [[0,5,1], [1,2,15], [2,4,7], [3,1,0], [5,0,0], [8,5,0], [13,2,0]],           #max, Hay cell 0
              [[0,5,0], [1,3,0], [2,5,1], [3,1,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],  #min, Hay cell 0
              [[0,5,1], [1,2,15], [2,4,7], [3,1,0], [5,0,0], [8,5,0], [13,3,0]],           #max, Hay cell 1
              [[0,5,0], [1,3,0], [2,5,1], [3,1,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],  #min, Hay cell 1
              [[0,5,1], [1,2,15], [2,4,7], [3,1,0], [5,0,0], [8,5,0], [13,2,0]],           #max, Hay cell 2
              [[0,5,0], [1,3,0], [2,5,1], [3,1,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],  #min, Hay cell 2
              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,5,0]],           #max, Hay cell 3
              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,4,0]],  #min, Hay cell 3
              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,1,0]],           #max, Hay cell 4
              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],  #min, Hay cell 4
              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,5,0]],           #max, Hay cell 5
              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,3,0]],  #min, Hay cell 5
              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,0,0]],           #max, Hay cell 6
              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]] ] #min, Hay cell 6                                                                                                                                
if len(sys.argv) > 1:
  icomb = int(float(sys.argv[1]))
combs = combs_all[icomb]

for icell in range(0,len(paramdicts)):
  print "icell = "+str(icell)
  if len(sys.argv) > 2 and int(float(sys.argv[2])) != icell:
    print "icell = "+str(icell)+", int(float(sys.argv[2])) = "+ str(int(float(sys.argv[2])))
    continue

  theseCoeffsAll = theseCoeffsAllAll[icell]

  spTimesThisVal = []
  spTimesThisVal2 = []
  ISIs_thismutval = []

  mytime = time.time()
  for iter in [0, 2, 5, 6, 8, -1]:
    if not exists('ifcurvesmut2_cs'+str(icell)+'_comb'+str(icomb)+'_iter'+str(iter)+'.sav'):
      print 'ifcurvesmut2_cs'+str(icell)+'_comb'+str(icomb)+'_iter'+str(iter)+'.sav does not exist, continuing'
      spTimesThisVal.append([])
      spTimesThisVal2.append([])
      ISIs_thismutval.append([])
      continue
    
    counter = -1

    nSpikes = []

    unpicklefile = open('ifcurvesmut2_cs'+str(icell)+'_comb'+str(icomb)+'_iter'+str(iter)+'.sav', 'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    spTimesThisVal.append(unpickledlist[1][0])
    spTimesThisVal2.append(unpickledlist[2][0])
    ISIs_thismutval.append(unpickledlist[3][0])

  picklelist = [ISIs_thismutval,spTimesThisVal,spTimesThisVal2,MT]
  file = open('ifcurvesmut2_cs'+str(icell)+'_comb'+str(icomb)+'.sav', 'w')
  pickle.dump(picklelist,file)
  file.close()

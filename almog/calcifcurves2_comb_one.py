#copied from calcifcurves2.py 30.10.2018
from neuron import h
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

h("""
load_file("myrun.hoc")
""")
print("printcounter="+str(printcounter)); printcounter = printcounter+1
h("""
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.001)

access a_soma
objref st1,syn1, sl
a_soma st1 = new IClamp(0.5)

double siteVec[2]
sl = new List()
sl=locateSites("apic",620)
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
apic[siteVec[0]] syn1 = new AlphaSynapse(siteVec[1])
//apic[41] syn1 = new AlphaSynapse(0.5)

syn1.onset = 3400
syn1.tau = 3
syn1.gmax = 0.0
syn1.e = 50

objref vsoma, vdend, tvec
vsoma = new Vector()
vdend = new Vector()
tvec = new Vector()
a_soma cvode.record(&v(0.5),vsoma,tvec)
apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)

v_init = -62
dt = 0.025
tstop = 1000
""")

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

print("printcounter="+str(printcounter)); printcounter = printcounter+1

for icell in range(0,len(paramdicts)):
  print "icell = "+str(icell)
  if len(sys.argv) > 2 and int(float(sys.argv[2])) != icell:
    print "icell = "+str(icell)+", int(float(sys.argv[2])) = "+ str(int(float(sys.argv[2])))
    continue

  paramdict = paramdicts[icell]
  print "Setting params..."
  setparams(paramdict)

  theseCoeffsAll = theseCoeffsAllAll[icell]

  spTimesThisVal = []
  spTimesThisVal2 = []
  ISIs_thismutval = []

  mytime = time.time()
  for iter in [0, 2, 5, 6, 8, -1]:
    if len(sys.argv) > 3 and int(float(sys.argv[3])) != iter:
      continue
    if exists('ifcurvesmut2_cs'+str(icell)+'_comb'+str(icomb)+'_iter'+str(iter)+'.sav'):
      print 'ifcurvesmut2_cs'+str(icell)+'_comb'+str(icomb)+'_iter'+str(iter)+'.sav exists, continuing'
      continue
    
    defVals = mutation_stuff.getdefvals()
    keyList = defVals.keys()
    counter = -1
    for igene in range(0,len(MT)):
     for imut in range(0,len(MT[igene])):
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
        isin = False
        for checkcomb in combs:
          if igene == checkcomb[0] and imut == checkcomb[1] and iallmutval == checkcomb[2]:
            isin = True
        if isin:
          mutval = allmutvals[iallmutval]
          spTimesThisVal = []
          spTimesThisVal2 = []
          ISIs_thismutval = []
  
          if iter >= 0:
            thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
          else:
            thisCoeff = 0
          mutText = ""
          for imutvar in range(0,len(MT[igene][imut])):
            if imutvar > 0 and imutvar%2==0:
              mutText = mutText+"\n"
            mutvars = allmutvars[iallmutval][imutvar]
            mutvals = allmutvals[iallmutval][imutvar]
            if type(mutvars) is str:
              mutvars = [mutvars]
            mutText = mutText + str(mutvars) + ": "
            for kmutvar in range(0,len(mutvars)):
              mutvar = mutvars[kmutvar]
              if (mutvars[kmutvar].find('off') > -1 and mutvars[kmutvar].find('offc') < 0) or mutvars[kmutvar].find('eh') > -1:
                newVal =  defVals[mutvar]+mutvals*thisCoeff
                defVals[mutvar] = defVals[mutvar]+mutvals*thisCoeff
                if mutvals >= 0 and kmutvar==0:
                  mutText = mutText + "+" + str(mutvals) +" mV"
                elif kmutvar==0:
                  mutText = mutText  + str(mutvals) +" mV"
              else:
                newVal = defVals[mutvar]*(mutvals**thisCoeff)
                defVals[mutvar] = defVals[mutvar]*(mutvals**thisCoeff)
                if kmutvar==0:
                  mutText = mutText + "*" + str(mutvals)
              if kmutvar < len(mutvars)-1:
                mutText = mutText + ", "
              mySuffix = mutvars[kmutvar][mutvars[kmutvar].find('_')+1:len(mutvars[kmutvar])]
              mySuffixInd = next((i for i,x in enumerate(mySuffixes) if x.find(mySuffix) > -1))
              isException = 0
              for jsuffe in range(0,len(mySuffixExceptions[mySuffixInd])):
                if mySuffixExceptions[mySuffixInd][jsuffe][0].find(mutvars[kmutvar]) > -1:
                  isException = 1
                  exceptionInd = jsuffe

              if not isException:
                print ("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mutvars[kmutvar]+""" = """+str(newVal))
                h("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mutvars[kmutvar]+""" = """+str(newVal))
              else:
                print ("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mySuffixExceptions[isuffix][j][1]+""" = """+str(newVal))
                h("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mySuffixExceptions[isuffix][j][1]+""" = """+str(newVal))

          print mutText

    thisCa = h.a_soma.cainf_cad

    spTimesThisCoeff = []
    spTimesThisCoeff2 = []
    ISIs = len(Is)*[0.0]
    nSpikes = []
    for iI in range(0,len(Is)):
          print "Running "+str(iI)
          tstop = 16000.0
          squareAmp = Is[iI]
          squareDur = 15800.0
          h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(thisCa)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 200
st1.dur = """+str(squareDur)+"""
syn1.gmax = 0
syn1.onset = 200 + """+str(BACdt)+""" 
  """)
          h.init()
          h.run()
  
          times=np.array(h.tvec)
          Vsoma=np.array(h.vsoma)
          Vdend=np.array(h.vdend)
          spikes = mytools.spike_times(times,Vsoma,-50,-50)
          spikes2 = mytools.spike_times(times,Vsoma,-50,inf)
          spTimesThisCoeff.append(spikes[:])
          spTimesThisCoeff2.append(spikes2[:])
          nSpikes1 = len(spikes)
          nSpikes2 = sum([1 for x in spikes if x >= 500.0])
          nSpikes.append(nSpikes2)
  
          if nSpikes1 > 5:
            spts = spikes[len(spikes)-5:len(spikes)]
            ISIs[iI] = mean([y-x for x,y in zip(spts[0:4],spts[1:5])])
          else:
            ISIs[iI] = 1.0e10
  
    defVals = mutation_stuff.getdefvals()
    keyList = defVals.keys()
    #Restore default values:
    for imutvar in range(0,len(MT[igene][imut])):
      mutvars = allmutvars[iallmutval][imutvar]
      mutvals = allmutvals[iallmutval][imutvar]
      if type(mutvars) is str:
        mutvars = [mutvars]
      for kmutvar in range(0,len(mutvars)):
        newVal = defVals[mutvars[kmutvar]]
        mySuffix = mutvars[kmutvar][mutvars[kmutvar].find('_')+1:len(mutvars[kmutvar])]
        mySuffixInd = next((i for i,x in enumerate(mySuffixes) if x.find(mySuffix) > -1))
        isException = 0
        for jsuffe in range(0,len(mySuffixExceptions[mySuffixInd])):
          if mySuffixExceptions[mySuffixInd][jsuffe][0].find(mutvars[kmutvar]) > -1:
            isException = 1
            exceptionInd = jsuffe
        if not isException:
          h("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mutvars[kmutvar]+""" = """+str(defVals[mutvars[kmutvar]]))
        else:
          h("""forall if(ismembrane(\""""+mySuffix+"""\")) """+mySuffixExceptions[isuffix][j][1]+""" = """+str(defVals[mutvars[kmutvar]]))
    spTimesThisVal.append(spTimesThisCoeff[:])
    spTimesThisVal2.append(spTimesThisCoeff2[:])
    ISIs_thismutval.append(ISIs[:])

    picklelist = [ISIs_thismutval,spTimesThisVal,spTimesThisVal2,MT]
    file = open('ifcurvesmut2_cs'+str(icell)+'_comb'+str(icomb)+'_iter'+str(iter)+'.sav', 'w')
    pickle.dump(picklelist,file)
    file.close()

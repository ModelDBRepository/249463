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

v0 = -62
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
mySuffixes = mutation_stuff.getsuffixes()
mySuffixExceptions = mutation_stuff.getsuffixexceptions()

unpicklefile = open('scalings_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
theseCoeffsAllAll = unpickledlist[0]
theseMutValsAllAll = unpickledlist[2]

paramdicts = []
paramdicts.append({'transvec.x(31)': 1.0, 'transvec.x(32)': 1.0, 'transvec.x(20)': 1.0, 'transvec.x(21)': 1.0, 'transvec.x(25)': 1.0, 'transvec.x(26)': 1.0}) # 4-6 spikes per burst, control
paramdicts.append({'transvec.x(31)': 1.25, 'transvec.x(32)': 1.25})                                                                                           # 4-5 spikes per burst         
paramdicts.append({'transvec.x(31)': 1.5, 'transvec.x(32)': 1.5})                                                                                             # 3-4 spikes per burst         
paramdicts.append({'transvec.x(31)': 2.0, 'transvec.x(32)': 2.0})                                                                                             # 3-4 spikes per burst        
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0})                                                                                             # 2-3 spikes per burst       
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0, 'transvec.x(20)': 1.3, 'transvec.x(21)': 1.3, 'transvec.x(25)': 1.3, 'transvec.x(26)': 1.3}) # 2 spikes per burst        
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0, 'transvec.x(20)': 1.6, 'transvec.x(21)': 1.6, 'transvec.x(25)': 1.6, 'transvec.x(26)': 1.6}) # 1-2 spikes per burst     

lensToStart = [450.0]


for istartdist in range(0,len(lensToStart)):
 startdist = lensToStart[istartdist]

 unpicklefile = open('synlocs'+str(startdist)+'.sav', 'r')
 unpickledlist = pickle.load(unpicklefile)
 unpicklefile.close()
 Nsyns = unpickledlist[0]
 synlocs = unpickledlist[3]
 startdist = int(startdist)

 gsAllAll = []

 for icell in range(0,7):
  gsAll = []
  theseCoeffsAll = theseCoeffsAllAll[icell]

  h("""
load_file("myrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.001)

access a_soma
objref st1,syn1, sl, syns["""+str(Nsyns)+"""]
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
""")
  paramdict = paramdicts[icell]
  setparams(paramdict)

  for istim in range(0,Nsyns):
    h("""
siteVec[0] = """+str(synlocs[istim][0])+"""
siteVec[1] = """+str(synlocs[istim][1])+"""
apic[siteVec[0]] {
  syns["""+str(istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(istim)+"""].e = 0
  syns["""+str(istim)+"""].tau = 5
  syns["""+str(istim)+"""].onset = 10000 + """+str(BACdt)+"""
}
""")

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
        if len(sys.argv) > 1 and int(float(sys.argv[1])) != counter:
          continue
        gsThisMutVal = []
        minNSpikesThisMutVal = []
        close("all")
        f, axarr = plt.subplots(2, 2)
        maxCac = 0
        maxCadc = 0
        if exists('thresholddistalamp'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav'):
          print 'thresholddistalamp'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav exists, continuing'
          continue
        for iter in [0, 2, 5, 6, 8, -1]:
          gsThisIter = []
          if iter >= 0:
            thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
          else:
            thisCoeff = 0
          if iter == -1 and (igene > 0 or imut > 0 or iallmutval > 0):
            continue # do the control only once!
          print "iter="+str(iter)+", thisCoeff="+str(thisCoeff)
          
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
              if (mutvar.find('off') > -1 and mutvar.find('offc') < 0) or mutvar.find('eh') > -1:
                newVal =  defVals[mutvar]+mutvals*thisCoeff
                if mutvals >= 0 and kmutvar==0:
                  mutText = mutText + "+" + str(mutvals) +" mV"
                elif kmutvar==0:
                  mutText = mutText  + str(mutvals) +" mV"
              else:
                newVal = defVals[mutvar]*(mutvals**thisCoeff)
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
          nextgs = [0.00,0.003,0.0015]
          hasSpiked = 0
          hasErred = 0
          minNSpikes = inf
          for iterg in range(0,ITERS+2):
              thisg = nextgs[min(iterg,2)]
              for istim in range(0,Nsyns):
                h("syns["+str(istim)+"].gmax = "+str(thisg))
              h("""
tstop = """+str(tstop)+"""
cai0_ca_ion = """+str(thisCa)+"""
v_init = """+str(v0)+"""
st1.amp = 0
st1.del = 0
st1.dur = 0
""")
              h.init()
              try:
                h.run()
              except RuntimeError:
                hasErred = 1
                print "Too large g!"
                if iterg == 1:
                  nextgs = [0.0,0.0015,0.00075]
                  continue
                else:
                  nextgs = [nextgs[0],nextgs[2],nextgs[0]+nextgs[2]]
                continue
  
              times=np.array(h.tvec)
              Vsoma=np.array(h.vsoma)
              spikes = mytools.spike_times(times,Vsoma,-50,-50)
              nSpikes1 = len(spikes)
              hasSpiked = hasSpiked or (nSpikes1 > 0)
              if nSpikes1 > 0 and nSpikes1 < minNSpikes:
                minNSpikes = nSpikes1
  
              print "iterg="+str(iterg)+" done, g="+str(thisg)+", "+str(nSpikes1)+" spikes"
              if iterg==0 and nSpikes1 > 0:
                print "Even zero g causes spiking!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)+", iter="+str(iter)+", spike at "+str(spikes[0])
                nextgs = [0.0,0.0,0.0]
                break
              if iterg==1 and not hasSpiked:
                print "No spiking with iterg==1, adding 900% to the current! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
                nextgs = [nextgs[0],10.0*nextgs[1],5*nextgs[min(iterg,2)]]
                continue

              if iterg>=2 and iterg < ITERS+2:
                if nSpikes1 > 0:
                  nextgs = [nextgs[0],nextgs[2],0.5*nextgs[0]+0.5*nextgs[2]]
                else:
                  nextgs = [nextgs[2],nextgs[1],0.5*nextgs[1]+0.5*nextgs[2]]

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
  
          gsThisMutVal.append(nextgs[2])
          minNSpikesThisMutVal.append(minNSpikes)
  
        gsThisMut.append(gsThisMutVal[:])
        picklelist = [theseCoeffsAll,gsThisMutVal,minNSpikesThisMutVal,MT]
        file = open('thresholddistalamp'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav', 'w')
        pickle.dump(picklelist,file)
        file.close()
      gsThisGene.append(gsThisMut[:])
    gsAll.append(gsThisGene[:])
  gsAllAll.append(gsAll[:])
  #picklelist = [theseCoeffsAll,gsThisAll,MT]
  #file = open('thresholddistalamp'+str(startdist)+'_'+str(counter)+'.sav', 'w')
  #pickle.dump(picklelist,file)
  #file.close()
  

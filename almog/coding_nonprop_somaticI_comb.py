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
import itertools

random.seed(1)

v0 = -62
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
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

lensToStart = [0.0 + x*200 for x in range(0,6)]
lensToEnd = [200.0 + x*200 for x in range(0,5)]+[1325.]
lenToStartBasal = 0.0
lenToEndBasal = 356.0
synconductance = 0.00002
if len(sys.argv) > 2:
  synconductance = float(sys.argv[2])
somaticI = 0.0
if len(sys.argv) > 4:
  somaticI = float(sys.argv[4])

for icell in range(0,7):
  gCoeffsAllAll = []
  synlocsAllAll = []
  Nsyns_all = []
  theseCoeffsAll = theseCoeffsAllAll[icell]
  for istartdist in range(0,len(lensToStart)):
    lenToStart = lensToStart[istartdist]
    lenToEnd = lensToEnd[istartdist]
    unpicklefile = open('synlocs'+str(lenToStart)+'-'+str(lenToEnd)+'.sav', 'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    Nsyns = unpickledlist[0]
    synlocsAll = unpickledlist[3]

    lenToStart = int(lenToStart)
    lenToEnd = int(lenToEnd)

    Nsyns_all.append(Nsyns)
    synlocsAllAll.append(synlocsAll[:])

  unpicklefile = open('basalsynlocs'+str(lenToStartBasal)+'-'+str(lenToEndBasal)+'.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  Nsyns = unpickledlist[0]
  synlocsAll = unpickledlist[3]
  Nsyns_all.append(Nsyns)
  synlocsAllAll.append(synlocsAll[:])


  coding_inputs = list(itertools.product([0, 1], repeat=7))

  if len(sys.argv) > 3 and int(float(sys.argv[3])) != icell:
    continue
  h("""
load_file("myrun.hoc")
objref cvode, sl
cvode = new CVode()
cvode.active(1)
cvode.atol(0.001)

access a_soma

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

objref vsoma, vdend, tvec, recSite, vdends, cadends, casoma, st1
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
tvec = new Vector()
cadends = new List()
vdends = new List()
objref ns, syns["""+str(7*Nsyns)+"""]

a_soma st1 = new IClamp(0.5)
a_soma cvode.record(&v(0.5),vsoma,tvec)
apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)

v_init = -62
st1.amp = """+str(somaticI)+"""
st1.del = 9900
st1.dur = 200
dt = 0.025
""")
  paramdict = paramdicts[icell]
  setparams(paramdict)

  for idist in range(0,6):
    h("""
cadends.append(new Vector())
vdends.append(new Vector())
sl = locateSites("apic","""+str((lensToStart[idist]+lensToEnd[idist])*0.5)+""")
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
apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadends.o[cadends.count()-1],tvec)
apic[siteVec[0]] cvode.record(&cai(siteVec[1]),vdends.o[vdends.count()-1],tvec)
""")
    for istim in range(0,Nsyns):
      h("""
siteVec[0] = """+str(synlocsAllAll[idist][istim][0])+"""
siteVec[1] = """+str(synlocsAllAll[idist][istim][1])+"""
apic[siteVec[0]] {
  syns["""+str(Nsyns*idist+istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(Nsyns*idist+istim)+"""].e = 0
  syns["""+str(Nsyns*idist+istim)+"""].tau = 5
  syns["""+str(Nsyns*idist+istim)+"""].onset = 10000
  syns["""+str(Nsyns*idist+istim)+"""].gmax = """+str(synconductance)+"""
}
""")
  h("""
cadends.append(new Vector())
vdends.append(new Vector())
sl = locateSites("dend","""+str((lenToStartBasal+lenToEndBasal)*0.5)+""")
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
dend[siteVec[0]] cvode.record(&cai(siteVec[1]),cadends.o[cadends.count()-1],tvec)
dend[siteVec[0]] cvode.record(&cai(siteVec[1]),vdends.o[vdends.count()-1],tvec)
""")
  for istim in range(0,Nsyns):
    h("""
siteVec[0] = """+str(synlocsAllAll[6][istim][0])+"""
siteVec[1] = """+str(synlocsAllAll[6][istim][1])+"""
dend[siteVec[0]] {     
  syns["""+str(Nsyns*6+istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(Nsyns*6+istim)+"""].e = 0 
  syns["""+str(Nsyns*6+istim)+"""].tau = 5
  syns["""+str(Nsyns*6+istim)+"""].onset = 10000
  syns["""+str(Nsyns*6+istim)+"""].gmax = """+str(synconductance)+"""
}
""")

  coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
  mytime = time.time()
  iitercounter = -1
  coding_outputs = []

  if exists('codings_nonprop'+str(synconductance)+'_cs'+str(icell)+'_comb'+str(icomb)+'_somaticI'+str(somaticI)+'.sav'):
    continue

  for iter in [0, 2, 6, 8, -1]:
    iitercounter = iitercounter + 1
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
          maxCac = 0
          maxCadc = 0

          if iter >= 0:
            thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
          else:
            thisCoeff = 0

          print "iter="+str(iter)+", thisCoeff="+str(thisCoeff)+", igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
          
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

    h("""
tstop = """+str(tstop)+"""
cai0_ca_ion = """+str(thisCa)+"""
v_init = """+str(v0)+"""
""")

    coding_outputs_thisiter = []
    for icoding in range(0,len(coding_inputs)):
      for idist in range(0,6):
        for istim in range(0,Nsyns):
          h("syns["+str(Nsyns*idist+istim)+"].gmax = "+str(synconductance*coding_inputs[icoding][idist]))
        print "syns["+str(Nsyns*idist)+"+istim].gmax = "+str(synconductance*coding_inputs[icoding][idist])
      for istim in range(0,Nsyns):
        h("syns["+str(Nsyns*6+istim)+"].gmax = "+str(synconductance*coding_inputs[icoding][6]))
      print "syns["+str(Nsyns*6)+"+istim].gmax = "+str(synconductance*coding_inputs[icoding][6])
      timenow = time.time()
      h.init()
      try:
        h.run()
      except RuntimeError:
        hasErred = 1
        print "Too large g!"

      times=np.array(h.tvec)
      Vsoma=np.array(h.vsoma)
      spikes = mytools.spike_times(times,Vsoma,-50,-50)
      nSpikes1 = len(spikes)
      coding_outputs_thisiter.append([max(array(h.cadends)[i]) for i in range(0,7)]+[nSpikes1]+[max(array(h.vdends)[i]) for i in range(0,7)])
    coding_outputs.append(coding_outputs_thisiter[:])

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


    picklelist = [[],coding_inputs,coding_outputs,[]]
    file = open('codings_nonprop'+str(synconductance)+'_cs'+str(icell)+'_comb'+str(icomb)+'_somaticI'+str(somaticI)+'.sav', 'w')
    pickle.dump(picklelist,file)
    file.close()


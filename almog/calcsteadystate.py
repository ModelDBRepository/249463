from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import time
import sys
from setparams import *
from os.path import exists

v0 = -62
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
fs = 8
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

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
paramdicts.append({})                                                                                                                                         # 4-6 spikes per burst, control
paramdicts.append({'transvec.x(31)': 1.25, 'transvec.x(32)': 1.25})                                                                                           # 4-5 spikes per burst         
paramdicts.append({'transvec.x(31)': 1.5, 'transvec.x(32)': 1.5})                                                                                             # 3-4 spikes per burst         
paramdicts.append({'transvec.x(31)': 2.0, 'transvec.x(32)': 2.0})                                                                                             # 3-4 spikes per burst        
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0})                                                                                             # 2-3 spikes per burst       
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0, 'transvec.x(20)': 1.3, 'transvec.x(21)': 1.3, 'transvec.x(25)': 1.3, 'transvec.x(26)': 1.3}) # 2 spikes per burst        
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0, 'transvec.x(20)': 1.6, 'transvec.x(21)': 1.6, 'transvec.x(25)': 1.6, 'transvec.x(26)': 1.6}) # 1-2 spikes per burst

spTimesAll = []
timesAll = []
VsomasAll = []
CasomasAll = []
VdendsAll = []
CadendsAll = []
ISIs_all = []

h("""
load_file("myrun.hoc")
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

objref vsoma, vdend, casoma, cadend, tvec
vsoma = new Vector()
vdend = new Vector()
casoma = new Vector()
cadend = new Vector()
tvec = new Vector()
a_soma cvode.record(&v(0.5),vsoma,tvec)
a_soma cvode.record(&cai(0.5),casoma,tvec)
apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)
apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend,tvec)

v_init = -62
dt = 0.025
tstop = 1000
""")                                          

for icell in range(0,7):
 paramdict = paramdicts[icell]
 setparams(paramdict)

 theseCoeffsAll = theseCoeffsAllAll[icell]
 spTimesAll = []
 spTimesAll2 = []
 ISIs_all = []

 counter = -1
 for igene in range(0,len(MT)):
  spTimesThisGene = []
  timesThisGene = []
  VsomasThisGene = []
  CasomasThisGene = []
  VdendsThisGene = []
  CadendsThisGene = []
  for imut in range(0,len(MT[igene])):
    spTimesThisMut = []
    timesThisMut = []
    VsomasThisMut = []
    CasomasThisMut = []
    VdendsThisMut = []
    CadendsThisMut = []
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
      spTimesThisMutVal = []
      timesThisMutVal = []
      VsomasThisMutVal = []
      CasomasThisMutVal = []
      VdendsThisMutVal = []
      CadendsThisMutVal = []
      for iter in [0, 2, 5, 6, 8, -1]:
        if iter >= 0:
          thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
        else:
          thisCoeff = 0
        if iter == -1 and (igene > 0 or imut > 0 or iallmutval > 0):
          continue # do the control only once!
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
        tstop = 4000.0
        squareAmp = 0.8
        squareDur = 3800.0
        h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 200
st1.dur = """+str(squareDur)+"""
syn1.gmax = 0
""")
        h.init()
        h.run()
  
        times=np.array(h.tvec)
        Vsoma=np.array(h.vsoma)
        Vdend=np.array(h.vdend)
        Casoma=np.array(h.casoma)
        Cadend=np.array(h.cadend)
        spikes = mytools.spike_times(times,Vsoma,-50,-50)
        spTimesThisCoeff = spikes[:]
        nSpikes1 = len(spikes)
  
        if nSpikes1 > 7:
          spts = spikes[len(spikes)-6:len(spikes)]
          istart = next((i for i,x in enumerate(times) if x > spts[0]))
          iend = next((i for i,x in enumerate(times) if x > spts[4]))+7
          nsteps = iend-istart-1
          tdiff = [y-x for x,y in zip(times[istart:iend-1],times[istart+1:iend])]
          cadiff = [y-x for x,y in zip(Casoma[istart:iend-1],Casoma[istart+1:iend])]
          caddiff = [y-x for x,y in zip(Cadend[istart:iend-1],Cadend[istart+1:iend])]
          caderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],cadiff[0:nsteps-1])]
          caderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],cadiff[1:nsteps])]
          caderiv = [(x+y)/2.0 for x,y in zip(caderiv1,caderiv2)]
          cadderiv = [y/x for x,y in zip(tdiff,caddiff)]
  
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

        spTimesThisMutVal.append(spTimesThisCoeff[:])
        timesThisMutVal.append(times[:])
        VsomasThisMutVal.append(Vsoma[:])
        CasomasThisMutVal.append(Casoma[:])
        VdendsThisMutVal.append(Vdend[:])
        CadendsThisMutVal.append(Cadend[:])
        if iter==-1:
          picklelist = [theseCoeffsAll,times,Vsoma,Casoma,Vdend,Cadend,MT]
          file = open('steadystate_cs'+str(icell)+'_control.sav', 'w')
          pickle.dump(picklelist,file)
          file.close()
  
      spTimesThisMut.append(spTimesThisMutVal[:])
      timesThisMut.append(timesThisMutVal[:])
      VsomasThisMut.append(VsomasThisMutVal[:])
      CasomasThisMut.append(CasomasThisMutVal[:])
      VdendsThisMut.append(VdendsThisMutVal[:])
      CadendsThisMut.append(CadendsThisMutVal[:])
      picklelist = [theseCoeffsAll,timesThisMutVal,VsomasThisMutVal,CasomasThisMutVal,VdendsThisMutVal,CadendsThisMutVal,MT]
      file = open('steadystate_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'w')
      pickle.dump(picklelist,file)
      file.close()
    spTimesThisGene.append(spTimesThisMut[:])
    timesThisGene.append(timesThisMut[:])
    VsomasThisGene.append(VsomasThisMut[:])
    CasomasThisGene.append(CasomasThisMut[:])
    VdendsThisGene.append(VdendsThisMut[:])
    CadendsThisGene.append(CadendsThisMut[:])
  
  spTimesAll.append(spTimesThisGene[:])
  timesAll.append(timesThisGene[:])
  VsomasAll.append(VsomasThisGene[:])
  CasomasAll.append(CasomasThisGene[:])
  VdendsAll.append(VdendsThisGene[:])
  CadendsAll.append(CadendsThisGene[:])
    
#picklelist = [theseCoeffsAll,VsomasAll,CasomasAll,VdendsAll,CadendsAll,MT]
#file = open('steadystate.sav', 'w')
#pickle.dump(picklelist,file)
#file.close()
  

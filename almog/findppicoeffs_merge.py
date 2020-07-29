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

#lensToStart = [100.0 + x*50 for x in range(0,16)]
lensToStart = [150.0, 300.0, 450.0, 600.0]

ISIs = [10.0*x for x in range(0,51)]
currCoeff = 1.1

gCoeffsAllAllAll = []

for istartdist in range(0,len(lensToStart)):
 startdist = lensToStart[istartdist]
 if len(sys.argv) > 2 and int(float(sys.argv[2])) != istartdist:
   continue

 unpicklefile = open('synlocs'+str(startdist)+'.sav', 'r')
 unpickledlist = pickle.load(unpicklefile)
 unpicklefile.close()
 Nsyns = unpickledlist[0]
 synlocs = unpickledlist[3]
 startdist = int(startdist)

 unpicklefile = open('thresholddistalamp'+str(startdist)+'.sav', 'r')
 unpickledlist = pickle.load(unpicklefile)
 unpicklefile.close()
 gsAllAll = unpickledlist[1]

 gCoeffsAllAll = []

 for icell in range(0,7):
  theseCoeffsAll = theseCoeffsAllAll[icell]
  gCoeffsAll = []
  if len(sys.argv) > 3 and int(float(sys.argv[3])) != icell:
    continue

  h("""
load_file("myrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.001)

access a_soma
objref st1,syn1, sl, syns["""+str(2*Nsyns)+"""]
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
  syns["""+str(istim)+"""].onset = 10000
  syns["""+str(Nsyns+istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(Nsyns+istim)+"""].e = 0
  syns["""+str(Nsyns+istim)+"""].tau = 5
  syns["""+str(Nsyns+istim)+"""].onset = 10000
}
""")

  coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

  counter = -1
  
  for igene in range(0,len(MT)):
    gCoeffsThisGene = []
    for imut in range(0,len(MT[igene])):
      gCoeffsThisMut = []
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
        gCoeffsThisMutVal = []
        close("all")
        f, axarr = plt.subplots(2, 2)
        maxCac = 0
        maxCadc = 0
        if exists('PPIcoeffs'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav'):
          print 'PPIcoeffs'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav exists, merging'
          unpicklefile = open('PPIcoeffs'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav', 'r')
          unpickledlist = pickle.load(unpicklefile)
          unpicklefile.close()
          #theseCoeffsAll = unpickledlist[0]
          gCoeffsThisMutVal = unpickledlist[1]
          #ISIs = unpickledlist[2]

        iters = [0, 2, 5, 6, 8, -1]
        nspsThisVal = []
        iitercounter = -1
        for iiter in range(0,len(iters)):          
          nspsThisIter = []
          iter = iters[iiter]
          gCoeffsThisIter = []

          if iter >= 0:
            thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
          else:
            thisCoeff = 0
          if iter == -1 and (igene > 0 or imut > 0 or iallmutval > 0):
            continue # do the control only once!
          if iter == 5:
            gCoeffsThisIter.append([])
            continue

          iitercounter = iitercounter + 1
          if len(gCoeffsThisMutVal) > iitercounter:
            print "counter="+str(counter)+", skipping iiter="+str(iiter)+" (iitercounter="+str(iitercounter)+") as len(gCoeffsThisMutVal)="+str(len(gCoeffsThisMutVal))
            continue

          thisg = gsAllAll[icell][igene][imut][iallmutval][iiter]
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
 
          for iISI in range(0,len(ISIs)):
            gCoeffsThisISI = []
            PPIdt = ISIs[iISI]
            nextCoeffs = [0,15.0,4.0]
            hasSpiked = 0
            for iterI in range(0,ITERS+2):
              for istim in range(0,Nsyns):
                h("syns["+str(istim)+"].gmax = "+str(thisg*currCoeff))
                h("syns["+str(Nsyns+istim)+"].gmax = "+str(thisg*currCoeff*nextCoeffs[min(iterI,2)]))
                h("syns["+str(Nsyns+istim)+"].onset = "+str(10000+PPIdt))
              h("""
tstop = """+str(tstop)+"""
cai0_ca_ion = """+str(thisCa)+"""
v_init = """+str(v0)+"""
st1.amp = 0
st1.del = 0
st1.dur = 0
""")
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
              print "nextCoeffs="+str(nextCoeffs)+", "+str(nSpikes1)+" spikes, simulation done in "+str(time.time()-timenow)+" seconds"
              #if iterI==0:
              #nSpikes_normal = 1
              nSpikes_normal = 1 + (icell<=2 and startdist<=400 or icell==3 and startdist<=300 or icell==4 and startdist<=250)
              hasSpiked = hasSpiked or (nSpikes1 > nSpikes_normal)
              if iterI == 0 and hasSpiked:
                print "istartdist="+str(istartdist)+", icell="+str(icell)+", igene="+str(igene)+", imut="+str(imut)+", iallmuval="+str(iallmutval)+", iiter="+str(iiter)+": extra spikes elicited for iterI=0!"
              if iterI > 0 and not hasSpiked:
                print "istartdist="+str(istartdist)+", icell="+str(icell)+", igene="+str(igene)+", imut="+str(imut)+", iallmuval="+str(iallmutval)+", iiter="+str(iiter)+": no extra spikes for iterI>0!"
                nextCoeffs = [nextCoeffs[1],2*nextCoeffs[1],1.5*nextCoeffs[1]]
                continue
              if iterI > 1 and nSpikes1 > nSpikes_normal:
                nextCoeffs = [nextCoeffs[0],nextCoeffs[2],0.5*(nextCoeffs[0]+nextCoeffs[2])]
              if iterI > 1 and nSpikes1 <= nSpikes_normal:
                nextCoeffs = [nextCoeffs[2],nextCoeffs[1],0.5*(nextCoeffs[2]+nextCoeffs[1])]              
              #print str(nSpikes1)+", nextCoeffs="+str(nextCoeffs)      
            gCoeffsThisISI = nextCoeffs[:]
            gCoeffsThisIter.append(gCoeffsThisISI[:])


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

          gCoeffsThisMutVal.append(gCoeffsThisIter[:])
          picklelist = [theseCoeffsAll,gCoeffsThisMutVal,ISIs,MT]
          file = open('PPIcoeffs'+str(startdist)+'_cs'+str(icell)+'_'+str(counter)+'.sav', 'w')
          pickle.dump(picklelist,file)
          file.close()
        gCoeffsThisMut.append(gCoeffsThisMutVal[:])
      gCoeffsThisGene.append(gCoeffsThisMut[:])
    gCoeffsAll.append(gCoeffsThisGene[:])
  gCoeffsAllAll.append(gCoeffsAll[:])
 gCoeffsAllAllAll.append(gCoeffsAllAll[:])
#picklelist = [theseCoeffsAll,gCoeffsThisAll,MT]
#file = open('thresholddistalamp'+str(startdist)+'_'+str(counter)+'.sav', 'w')
#pickle.dump(picklelist,file)
#file.close()
  

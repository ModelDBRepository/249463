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

import mutation_stuff
MT = mutation_stuff.getMT()
defVals = mutation_stuff.getdefvals()
keyList = defVals.keys()
mySuffixes = mutation_stuff.getsuffixes()
mySuffixExceptions = mutation_stuff.getsuffixexceptions()

indIsLc = 2
nSpikesPerLc = 5

theseCoeffsAllAll = []

##squareAmps = [0.856,0,0.802]
##epsp_gmaxs = [0,0.0392,0.0367]
#squareAmps = [0.781,0,0.699]
#epsp_gmaxs = [0,0.0376,0.0337]

BACdt = 2.5

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

syn1.onset = 3400+"""+str(BACdt)+"""
syn1.tau = 3
syn1.gmax = 0.005
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

ITERS = 30

paramdicts = []
paramdicts.append({})                                                                                                                                         # 4-6 spikes per burst, control
paramdicts.append({'transvec.x(31)': 1.25, 'transvec.x(32)': 1.25})                                                                                           # 4-5 spikes per burst
paramdicts.append({'transvec.x(31)': 1.5, 'transvec.x(32)': 1.5})                                                                                             # 3-4 spikes per burst
paramdicts.append({'transvec.x(31)': 2.0, 'transvec.x(32)': 2.0})                                                                                             # 3-4 spikes per burst
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0})                                                                                             # 2-3 spikes per burst
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0, 'transvec.x(20)': 1.3, 'transvec.x(21)': 1.3, 'transvec.x(25)': 1.3, 'transvec.x(26)': 1.3}) # 2 spikes per burst  
paramdicts.append({'transvec.x(31)': 4.0, 'transvec.x(32)': 4.0, 'transvec.x(20)': 1.6, 'transvec.x(21)': 1.6, 'transvec.x(25)': 1.6, 'transvec.x(26)': 1.6}) # 1-2 spikes per burst

theseCoeffsAllAll = []
theseMutValsAllAll = []
theseMutVarsAllAll = []

for icell in range(0,len(paramdicts)):

 paramdict = paramdicts[icell]
 setparams(paramdict)

 theseCoeffsAll = []
 theseMutValsAll = []
 theseMutVarsAll = []

 unpicklefile = open('controlamps_cs'+str(icell)+'.sav', 'r')
 unpickledlist = pickle.load(unpicklefile)
 unpicklefile.close()

 mySquareAmpsAll = unpickledlist[0]
 mySynAmpsAll = unpickledlist[1]

 squareAmps = [mySquareAmpsAll[i][0][-1] for i in range(0,len(mySquareAmpsAll))]
 synAmps = [mySynAmpsAll[i][0][-1] for i in range(0,len(mySynAmpsAll))]
 squareDurs = [5,5,5]


 unpicklefile = open('control_cs'+str(icell)+'.sav', 'r')
 unpickledlist = pickle.load(unpicklefile)
 unpicklefile.close()

 Is = [0.7,0.75,0.8,0.85]
 Is_control = unpickledlist[19]
 spikfreqs_control = mytools.interpolate(Is_control,unpickledlist[0],Is)
 timesc_control = unpickledlist[1] 
 Vsomac_control = unpickledlist[2]
 VDerivc_control = unpickledlist[3]
 VDcoeff_control = unpickledlist[4][0]



 counter = -1
 for igene in range(0,len(MT)):
  #for igene in range(0,1):
  theseCoeffsGene = []
  for imut in range(0,len(MT[igene])):
    #for imut in range(0,1):
    theseCoeffsMut = []
    nVals = len(MT[igene][imut])*[0]
    thesemutvars = []
    for imutvar in range(0,len(MT[igene][imut])):
      thesemutvars.append(MT[igene][imut][imutvar][0])
      if type(MT[igene][imut][imutvar][1]) is int or type(MT[igene][imut][imutvar][1]) is float:
        MT[igene][imut][imutvar][1] = [MT[igene][imut][imutvar][1]]
      nVals[imutvar] = len(MT[igene][imut][imutvar][1])
    cumprodnVals = cumprod(nVals)
    allmutvars = cumprodnVals[len(MT[igene][imut])-1]*[thesemutvars[:]]
    allmutvals = []
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      allmutvals.append([0]*len(thesemutvars))
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      for imutvar in range(0,len(MT[igene][imut])):
        if imutvar==0:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][iallmutval%nVals[imutvar]]
        else:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][(iallmutval/cumprodnVals[imutvar-1])%nVals[imutvar]]
    theseMutValsAll.append(allmutvals[:])  
    theseMutVarsAll.append(allmutvars[:])  
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      counter = counter + 1
      if len(sys.argv) > 1 and int(float(sys.argv[1])) != counter:
        continue
      if exists('scalings_cs'+str(icell)+'_'+str(counter)+'.sav'):
        continue
      nextCoeffs = [0.0,2.0,1.0]
      for iter in range(0,ITERS+2+3):
        thisCoeff = nextCoeffs[min(iter,2)]
   
        mutText = ""
        for imutvar in range(0,len(MT[igene][imut])):
          if imutvar > 0 and imutvar%2==0:
            mutText = mutText+"\n"
          mutvars = allmutvars[iallmutval][imutvar]
          if type(mutvars) is str:
            mutvars = [mutvars]
          mutText = mutText + str(mutvars) + ": "
          mutvals = allmutvals[iallmutval][imutvar]
          #if type(mutvals) is list:
          #  mutvals = max([mutvals[max(range(len(mutvals)), key=lambda i: mutvals[i])],-mutvals[min(range(len(mutvals)), key=lambda i: mutvals[i])]])
          for kmutvar in range(0,len(mutvars)):
            if (mutvars[kmutvar].find('off') > -1 and mutvars[kmutvar].find('offc') < 0) or mutvars[kmutvar].find('eh') > -1:
              newVal = defVals[mutvars[kmutvar]]+thisCoeff*mutvals
              if mutvals >= 0 and kmutvar==0:
                mutText = mutText + "+" + str(mutvals*thisCoeff) +" mV"
              elif kmutvar==0:
                mutText = mutText  + str(mutvals*thisCoeff) +" mV"
            else:
              newVal = defVals[mutvars[kmutvar]]*(mutvals**thisCoeff)
              if kmutvar==0:
                mutText = mutText + "*" + str(mutvals**thisCoeff)
            if kmutvar < len(mutvars)-1:
              mutText = mutText + ", "
            #else:
            #  mutText = mutText + "\n"

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
    
        close("all")

        f, axarr = plt.subplots(2, 4)
        for ix in range(0,3):
          for iy in range(0,2):
            axarr[iy,ix].set_position([0.05+0.3*ix, 0.05+0.4*(1-iy), 0.23, 0.3])

        tstop = 4000.0
        ampCoeffs = [0.85,1.15]
        isChanged = False
        for icond in range(0,3):
          ############################################# Condition 1-3: Response to short stimuli #############################################
          for iampCoeff in range(0,len(ampCoeffs)):
            squareAmp = squareAmps[icond]*ampCoeffs[iampCoeff]
            squareDur = squareDurs[icond]
            epsp_gmax = synAmps[icond]*ampCoeffs[iampCoeff]
            h("""
tstop = """+str(tstop)+"""
v_init = -62
cai0_ca_ion = """+str(thisCa)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 3400
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
""")
            h.init()
            h.run()
            times=np.array(h.tvec)
            Vdend=np.array(h.vdend)
            Vsoma=np.array(h.vsoma)
            spikes = mytools.spike_times(times,Vsoma,-20,-45)
            nSpikes1 = len(spikes)
            print "icond="+str(icond)+",iamp="+str(iampCoeff)+", nSpikes1="+str(nSpikes1)+", st1.amp="+str(h.st1.amp)+", syn1.gmax="+str(h.syn1.gmax)

            axarr[iampCoeff,icond].plot(times, Vsoma)
            axarr[iampCoeff,icond].set_title("nspikes="+str(nSpikes1))
            axarr[iampCoeff,icond].set_xlim([3390,3500])
            axarr[iampCoeff,icond].set_ylim([-100,40])
            f.savefig("testtesttest.png")

            isChanged = isChanged or nSpikes1 > 0 and iampCoeff == 0 or nSpikes1 == 0 and iampCoeff == 1

        print isChanged
        if isChanged:
          if iter==0:
            print "Even null mutation causes different spiking!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
            nextCoeffs = [nextCoeffs[0],nextCoeffs[1],nextCoeffs[0]]
            f.suptitle(mutText)
            f.savefig("vrecs_mut"+str(igene)+"_"+str(imut)+"_"+str(iallmutval)+"_ITER"+str(iter)+".png")
            break
          if iter>=2 and iter < ITERS+2:
            nextCoeffs = [nextCoeffs[0],nextCoeffs[2],0.5*nextCoeffs[0]+0.5*nextCoeffs[2]]
          f.savefig("vrecs_mut"+str(igene)+"_"+str(imut)+"_"+str(iallmutval)+"_ITER"+str(iter)+".png")
          continue

        ############################################# Condition 4: IF curve #############################################
        spikfreqs = len(Is)*[0]
        for iI in range(0,len(Is)):
          tstop = 7200.0
          squareAmp = Is[iI]
          squareDur = 3800.0
          epsp_gmax = 0.0
          h("""
tstop = """+str(tstop)+"""
v_init = -62
cai0_ca_ion = """+str(thisCa)+"""
st1.amp = """+str(squareAmp)+"""
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
""")
          h.init()
          h.run()

          times=np.array(h.tvec)
          Vsoma=np.array(h.vsoma)
          spikes = mytools.spike_times(times,Vsoma,-20,-45)
          spikfreqs[iI] = sum([1 for x in spikes if x >= 3700.0])/3.5
          if iI==2: # use the memb. pot. time course of 0.8nA for the limit cycle
            times_lc = times[:]
            Vsoma_lc = Vsoma[:]
            spikes_lc = spikes[:]
          f2, axarr2 = subplots(1, 1)
          axarr2.plot(times, Vsoma)
          axarr2.set_title("Perisomatic firing, nspikes="+str(len(spikes))+", after 300ms nspikes="+str(sum([1 for x in spikes if x >= 3700.0])))
          axarr2.set_ylim([-100,40])
          axarr2.set_xlim([3390,7200])
          f2.savefig("testscaling_"+str(iI)+".eps")

        spikfreqdiffsum = sum([abs(x-y) for x,y in zip(spikfreqs,spikfreqs_control)])
        spikfreqdiffrel = spikfreqdiffsum/sum(spikfreqs_control)
        axarr[0,3].plot(Is, spikfreqs)
        axarr[0,3].set_title("IF, diff="+str(spikfreqdiffrel))
        axarr[0,3].set_xlim([0,1.25])
        axarr[0,3].set_ylim([0,20])

        ############################################# Condition 5: Limit cycles ############################################

        print "spikfreqdiffrel="+str(spikfreqdiffrel)

        f.suptitle(mutText)
        if iter < ITERS+2:
          f.savefig("vrecs_mut"+str(igene)+"_"+str(imut)+"_"+str(iallmutval)+"_ITER"+str(iter)+".png")
        else:
          f.savefig("vrecs_mut"+str(igene)+"_"+str(imut)+"_"+str(iallmutval)+"_TEST"+str(iter-ITERS-2)+".png")

        isChanged = isChanged or spikfreqdiffrel > 0.1
        print isChanged
        if iter==0 and isChanged:
          print "Even null mutation causes different spiking!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
          nextCoeffs = [nextCoeffs[0],nextCoeffs[1],nextCoeffs[0]]
          break
        if iter==1 and not isChanged:
          print "This mutation effect does not alter spiking even when doubled!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
          continue
        if iter>=2 and iter < ITERS+2:
          if isChanged:
            nextCoeffs = [nextCoeffs[0],nextCoeffs[2],0.5*nextCoeffs[0]+0.5*nextCoeffs[2]]
          else:
            nextCoeffs = [nextCoeffs[2],nextCoeffs[1],0.5*nextCoeffs[1]+0.5*nextCoeffs[2]]
        if iter == ITERS+1:
          nextCoeffs = [nextCoeffs[2],nextCoeffs[2],nextCoeffs[2]*0.99]
        if iter == ITERS+2:
          nextCoeffs = [nextCoeffs[0],nextCoeffs[0],nextCoeffs[0]*1.0]
        if iter == ITERS+3:
          nextCoeffs = [nextCoeffs[0],nextCoeffs[0],nextCoeffs[0]*1.01]
      

      #Restore default values:
      for imutvar in range(0,len(MT[igene][imut])):
        mutvars = allmutvars[iallmutval][imutvar]
        if type(mutvars) is str:
          mutvars = [mutvars]
        mutvals = allmutvals[iallmutval][imutvar]
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


      theseCoeffsMut.append(nextCoeffs[0]+0.0)
      picklelist = [nextCoeffs[0]+0.0,igene,imut,iallmutval,counter,MT]
      file = open('scalings_cs'+str(icell)+'_'+str(counter)+'.sav', 'w')
      pickle.dump(picklelist,file)
      file.close()

    theseCoeffsGene.append(theseCoeffsMut[:])
  theseCoeffsAll.append(theseCoeffsGene[:])
 theseCoeffsAllAll.append(theseCoeffsAll[:])

#picklelist = [theseCoeffsAllAll,theseMutVarsAll,theseMutValsAll,MT]
#file = open('scalings.sav', 'w')
#pickle.dump(picklelist,file)
#file.close()

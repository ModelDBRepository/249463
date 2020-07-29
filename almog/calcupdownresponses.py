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

v0 = -80
ca0 = 0.0001
proximalpoints = [400]#,200]
distalpoints =   [850]#,900]
BACdt = 5.0
fs = 8
tstop = 5000.0
#epspdts = [0.25*x for x in range(-80,81)]
epspdts = [4.0*x for x in range(-40,-20)]+[2.0*x for x in range(-40,-20)]+[1.0*x for x in range(-40,-20)]+[0.5*x for x in range(-40,41)]+[1.0*x for x in range(21,41)]+[2.0*x for x in range(21,41)]+[4.0*x for x in range(21,41)]
epspdts_savetimecourses = [-40,-20,-15,-10,-5,0,5,10,15,20,40]

Is_st2 = 0.462
st2coeff = 0.40       #Somatic 5ms pulse
st2coeff_down = 1.35  #Somatic 5ms pulse
st1coeff = 0.85       #Proximal apical 200ms pulse
syn1coeff = 0.15      #Synaptic epsp-like input
syn1coeff_down = 0.15 #Synaptic epsp-like input

unpicklefile = open('apicalthresholds.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
IsAllAll_st1 = unpickledlist[0]
dists = unpickledlist[1]

unpicklefile = open('apicalthresholds_epsp.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
IsAllAll_syn1 = unpickledlist[0]

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

VsomaupAllAll = []
VsomadownAllAll = []
VdendupAllAll = []
VdenddownAllAll = []

counter = -1
for icell in range(0,7):
  theseCoeffsAll = theseCoeffsAllAll[icell]

  for idist in range(0,len(proximalpoints)):
    counter = counter + 1
    if len(sys.argv) > 2 and int(float(sys.argv[2])) != counter:
      continue
    
    proximalpoint = proximalpoints[idist]
    distalpoint = distalpoints[idist]
    fixedpoint = 700
    idist_proximal = dists.index(proximalpoint)
    idist_distal = dists.index(distalpoint)
    Is_syn1 = IsAllAll_syn1[icell][idist_distal]
    if type(Is_syn1) is list:
      Is_syn1 = Is_syn1[-1]
    Is_st1 = IsAllAll_st1[icell][idist_proximal]
    if type(Is_st1) is list:
      Is_st1 = Is_st1[-1]
    h("""
load_file("myrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.001)

access a_soma
objref st1,st2,syn1, sl
objref vsoma, casoma, sksoma, vdend, vdend2, vdend3, cadend3, skdend3, tvec, isyn
tvec = new Vector()
vsoma = new Vector()
casoma = new Vector()
sksoma = new Vector()
vdend = new Vector()
vdend2 = new Vector()
vdend3 = new Vector()
cadend3 = new Vector()
skdend3 = new Vector()
isyn = new Vector()

a_soma st2 = new IClamp(0.5)
st2.amp = """+str(st2coeff*Is_st2)+"""
st2.del = 1000
st2.dur = 5

double siteVec[2]
sl = new List()
sl=locateSites("apic","""+str(distalpoint)+""")
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
apic[siteVec[0]] syn1 = new epsp(siteVec[1])
//apic[41] syn1 = new AlphaSynapse(0.5)
syn1.tau0 = 0.5
syn1.tau1 = 5
syn1.onset = 1000 + """+str(BACdt)+"""
syn1.imax = """+str(syn1coeff*Is_syn1)+"""
apic[siteVec[0]] cvode.record(&syn1.i,isyn,tvec)
apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)

sl=locateSites("apic","""+str(proximalpoint)+""")
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
apic[siteVec[0]] st1 = new IClamp(siteVec[1])
st1.amp = """+str(st1coeff*Is_st1)+"""
st1.del = 700
st1.dur = 600
apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend2,tvec)

sl=locateSites("apic","""+str(fixedpoint)+""")
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
apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend3,tvec)
apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend3,tvec)
apic[siteVec[0]] cvode.record(&ik_sk(siteVec[1]),skdend3,tvec)
a_soma cvode.record(&v(0.5),vsoma,tvec)
a_soma cvode.record(&cai(0.5),casoma,tvec)
a_soma cvode.record(&ik_sk(0.5),sksoma,tvec)

v_init = -62
dt = 0.025
""")
    paramdict = paramdicts[icell]
    setparams(paramdict)

    styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
    #cols = ['#00aaaa','#11cc44','#55ee00','#bbaa00','#ee6600','#ff0000', '#aa00aa','#772277','#333333']
    cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
    coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
    lw = 0.5

    VsomaupAll = []
    VsomadownAll = []
    VdendupAll = []
    VdenddownAll = []

    mutcounter = -1

    for igene in range(0,len(MT)):
     VsomaupThisGene = []
     VsomadownThisGene = []
     VdendupThisGene = []
     VdenddownThisGene = []

     for imut in range(0,len(MT[igene])): 
      VsomaupThisMut = []
      VsomadownThisMut = []
      VdendupThisMut = []
      VdenddownThisMut = []
      spikesupThisMut = []
      spikesdownThisMut = []

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
        mutcounter = mutcounter + 1                                                                                                                                                               
        if len(sys.argv) > 1 and int(float(sys.argv[1])) != mutcounter:
          continue
        if exists('updownresponse_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav'):
          print 'updownresponse_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav exists'
          continue
        CasomaupThisMutVal = []
        CasomadownThisMutVal = []
        SKsomaupThisMutVal = []
        SKsomadownThisMutVal = []
        VdendupThisMutVal = []
        VdenddownThisMutVal = []
        Vdend2upThisMutVal = []
        Vdend2downThisMutVal = []
        Vdend3upThisMutVal = []
        Vdend3downThisMutVal = []
        Cadend3upThisMutVal = []
        Cadend3downThisMutVal = []
        SKdend3upThisMutVal = []
        SKdend3downThisMutVal = []
        spikesupThisMutVal = []
        spikesdownThisMutVal = []
        timecoursedataThisMutVal = []
        iters = [0,2,5,6,8,-1]
        #iters = [-1]
        for iiter in range(0,len(iters)):
          iter = iters[iiter]
          CasomaupThisIter = []
          CasomadownThisIter = []
          SKsomaupThisIter = []
          SKsomadownThisIter = []
          VdendupThisIter = []
          VdenddownThisIter = []
          Vdend2upThisIter = []
          Vdend2downThisIter = []
          Vdend3upThisIter = []
          Vdend3downThisIter = []
          Cadend3upThisIter = []
          Cadend3downThisIter = []
          SKdend3upThisIter = []
          SKdend3downThisIter = []
          spikesupThisIter = []
          spikesdownThisIter = []
          timecoursedataThisIter = []

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
          thisCa = h.a_soma.cainf_cad

          for idt in range(0,len(epspdts)):
            for iup in range(0,2):
              print "st1.amp = "+str(st1coeff*Is_st1*(1-iup))
              print "st2.amp = "+str((st2coeff*(1-iup) + st2coeff_down*iup)*Is_st2)
              print "syn1.imax = "+str((syn1coeff*(1-iup) + syn1coeff_down*iup)*Is_syn1)
              h("st1.amp = "+str(st1coeff*Is_st1*(1-iup)))
              h("st2.amp = "+str((st2coeff*(1-iup) + st2coeff_down*iup)*Is_st2))
              h("syn1.imax = "+str((syn1coeff*(1-iup) + syn1coeff_down*iup)*Is_syn1))
              h("""
tstop = """+str(tstop)+"""
cai0_ca_ion = """+str(thisCa)+"""
v_init = """+str(v0)+"""
syn1.onset = """+str(1000+epspdts[idt])+"""
""")
              h.init()
              try:
                h.run()
              except RuntimeError:
                print "Too large I!"
                continue
              times=np.array(h.tvec)
              Vsoma=np.array(h.vsoma)
              Casoma=np.array(h.casoma)
              SKsoma=np.array(h.sksoma)
              Vdend=np.array(h.vdend)
              Vdend2=np.array(h.vdend2)
              Vdend3=np.array(h.vdend3)
              Cadend3=np.array(h.cadend3)
              SKdend3=np.array(h.skdend3)
              spikes = mytools.spike_times(times,Vsoma,-35,-45)

              if iup == 1:
                CasomadownThisIter.append(max(Casoma))
                SKsomadownThisIter.append(max(SKsoma))
                VdenddownThisIter.append(max(Vdend))
                Vdend2downThisIter.append(max(Vdend2))
                Vdend3downThisIter.append(max(Vdend3))
                Cadend3downThisIter.append(max(Cadend3))
                SKdend3downThisIter.append(max(SKdend3))
                spikesdownThisIter.append(spikes[:])
              else:
                CasomaupThisIter.append(max(Casoma))
                SKsomaupThisIter.append(max(SKsoma))
                VdendupThisIter.append(max(Vdend))
                Vdend2upThisIter.append(max(Vdend2))
                Vdend3upThisIter.append(max(Vdend3))
                Cadend3upThisIter.append(max(Cadend3))
                SKdend3upThisIter.append(max(SKdend3))
                spikesupThisIter.append(spikes[:])
              if epspdts[idt] in epspdts_savetimecourses:
                times_tc = [680+x for x in range(0,641)]
                Vsoma_tc = mytools.interpolate(times,Vsoma,times_tc) 
                Casoma_tc = mytools.interpolate(times,Casoma,times_tc) 
                SKsoma_tc = mytools.interpolate(times,SKsoma,times_tc) 
                Vdend_tc = mytools.interpolate(times,Vdend,times_tc) 
                Vdend2_tc = mytools.interpolate(times,Vdend2,times_tc) 
                Vdend3_tc = mytools.interpolate(times,Vdend3,times_tc) 
                Cadend3_tc = mytools.interpolate(times,Cadend3,times_tc) 
                SKdend3_tc = mytools.interpolate(times,SKdend3,times_tc) 
                picklelist = [times_tc[:],Vsoma_tc[:],Casoma_tc[:],SKsoma_tc[:],Vdend_tc[:],Vdend2_tc[:],Vdend3_tc[:],Cadend3_tc[:],SKdend3_tc[:],epspdts_savetimecourses]
                timecoursedataThisIter.append(picklelist[:])
                file = open('updownresponsetimecourse_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'w')
                pickle.dump(picklelist,file)
                file.close()
          
              
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

          CasomaupThisMutVal.append(CasomaupThisIter)
          SKsomaupThisMutVal.append(SKsomaupThisIter)
          VdendupThisMutVal.append(VdendupThisIter)
          Vdend2upThisMutVal.append(Vdend2upThisIter)
          Vdend3upThisMutVal.append(Vdend3upThisIter)
          Cadend3upThisMutVal.append(Cadend3upThisIter)
          SKdend3upThisMutVal.append(SKdend3upThisIter)
          CasomadownThisMutVal.append(CasomadownThisIter)
          SKsomadownThisMutVal.append(SKsomadownThisIter)
          VdenddownThisMutVal.append(VdenddownThisIter)
          Vdend2downThisMutVal.append(Vdend2downThisIter)
          Vdend3downThisMutVal.append(Vdend3downThisIter)
          Cadend3downThisMutVal.append(Cadend3downThisIter)
          SKdend3downThisMutVal.append(SKdend3downThisIter)
          spikesupThisMutVal.append(spikesupThisIter)
          spikesdownThisMutVal.append(spikesdownThisIter)
          timecoursedataThisMutVal.append(timecoursedataThisIter[:])
          picklelist = [theseCoeffsAllAll,timecoursedataThisMutVal,epspdts_savetimecourses,MT]
          file = open('updownresponsetimecourse_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'w')
          pickle.dump(picklelist,file)
          file.close()

        picklelist = [theseCoeffsAllAll,CasomaupThisMutVal,SKsomaupThisMutVal,VdendupThisMutVal,Vdend2upThisMutVal,Vdend3upThisMutVal,Cadend3upThisMutVal,SKdend3upThisMutVal,
                      CasomadownThisMutVal,SKsomadownThisMutVal,VdenddownThisMutVal,Vdend2downThisMutVal,Vdend3downThisMutVal,Cadend3downThisMutVal,SKdend3downThisMutVal,
                      spikesupThisMutVal,spikesdownThisMutVal,epspdts,MT]
        file = open('updownresponse_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'w')
        pickle.dump(picklelist,file)
        file.close()

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

morphology_file = "morphologies/cell1.asc"
biophys_file = "models/L5PCbiophys3.hoc"
template_file = "models/L5PCtemplate.hoc"
v0 = -80
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
for idefval in range(0,len(keyList)):
  if type(defVals[keyList[idefval]]) is not list:
    defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments
unpicklefile = open('scalings_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
theseCoeffsAllAll = unpickledlist[0]
theseMutValsAll = unpickledlist[2]

paramdicts = []
paramdicts.append({})                                                         # 1 spike per burst, control
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 1.6})                               # 1-2 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2})                               # 2-3 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.9})   # 3-4 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.625}) # 3-5 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.5})   # 4-6 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.3})   # 5-9 spikes per burst      

spTimesAllAll = []
timesAllAll = []
VsomasAllAll = []
CasomasAllAll = []
VdendsAllAll = []
CadendsAllAll = []

for icell in range(0,7):
  theseCoeffsAll = theseCoeffsAllAll[icell]
  spTimesAll = []
  spTimesAll2 = []
  ISIs_all = []
  morphology_file = "morphologies/cell1.asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
objref st1
st1 = new IClamp(0.5)
L5PC.soma st1
forsec L5PC.somatic {
}
forsec L5PC.apical {
}
L5PC.distribute_channels("apic","gIhbar_Ih",2,-0.8696,3.6161,0.0,1.0*2.0870,0.0002)
L5PC.distribute_channels("apic","gCa_HVAbar_Ca_HVA",3,1.0,0.1,685.0,885.0,1.0*0.000555)
L5PC.distribute_channels("apic","gCa_LVAstbar_Ca_LVAst",3,1.0,0.01,685.0,885.0,1.0*0.0187)
objref sl,st2,ns,syn1,con1,isyn, tvec
isyn = new Vector()
tvec = new Vector()
sl = new List()
double siteVec[2]
sl = L5PC.locateSites("apic","""+str(distalpoint)+""")
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
print "distalpoint gCa_HVA: ", L5PC.apic[siteVec[0]].gCa_HVAbar_Ca_HVA
print "distalpoint gCa_LVA: ", L5PC.apic[siteVec[0]].gCa_LVAstbar_Ca_LVAst
access L5PC.apic[siteVec[0]]
st2 = new IClamp(siteVec[1])
st2.amp = 0
L5PC.apic[siteVec[0]] {
  st2
  syn1 = new epsp(siteVec[1])
  syn1.tau0 = 0.5
  syn1.tau1 = 5
  syn1.onset = 145 + """+str(BACdt)+""" 
  cvode.record(&syn1.i,isyn,tvec)
}
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
access L5PC.soma
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
access L5PC.apic[siteVec[0]]
cvode.record(&v(siteVec[1]),vdend,tvec)
cvode.record(&cai(siteVec[1]),cadend,tvec)
sl = new List()
sl = L5PC.locateSites("apic","""+str(proximalpoint)+""")
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
print "proximalpoint gCa_HVA: ", L5PC.apic[siteVec[0]].gCa_HVAbar_Ca_HVA
print "proximalpoint gCa_LVA: ", L5PC.apic[siteVec[0]].gCa_LVAstbar_Ca_LVAst
access L5PC.apic[siteVec[0]]
recSite = new IClamp(siteVec[1])
recSite.amp = 0
L5PC.apic[siteVec[0]] {
        recSite
}
access L5PC.apic[siteVec[0]]
cvode.record(&v(siteVec[1]),vdend2,tvec)
cvode.record(&cai(siteVec[1]),cadend2,tvec)
access L5PC.soma
isoma = new Vector()
cvode.record(&st1.i,isoma,tvec)
""")

  paramdict = paramdicts[icell]
  setparams(paramdict)  

  spTimesAll = []
  timesAll = []
  VsomasAll = []
  CasomasAll = []
  VdendsAll = []
  CadendsAll = []
  
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
            if mutvar.find('offm') > -1 or mutvar.find('offh') > -1 or mutvar.find('ehcn') > -1:
              newVal =  [x+mutvals*thisCoeff for x in defVals[mutvar]]
              if mutvals >= 0 and kmutvar==0:
                mutText = mutText + "+" + str(mutvals) +" mV"
              elif kmutvar==0:
                mutText = mutText  + str(mutvals) +" mV"
            else:
              newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
              if kmutvar==0:
                mutText = mutText + "*" + str(mutvals)
            if kmutvar < len(mutvars)-1:
              mutText = mutText + ", "
            if mutvar.find('_Ih') > -1:
              updateThese = [1,1,1]
            elif mutvar.find('_Ca_HVA') > -1 or mutvar.find('_Ca_LVAst') > -1 or mutvar.find('_SKv3.1') > -1 or mutvar.find('_Ca_HVA') > -1 or mutvar.find('_SK_E2') > -1 or mutvar.find('_NaTa_t') > -1 or mutvar.find('_CaDynamics_E2') > -1:
              updateThese = [1,1,0]
            elif mutvar.find('_K_Pst') > -1 or mutvar.find('_K_Tst') > -1 or mutvar.find('_Nap_Et2') > -1:
              updateThese = [1,0,0]
            elif mutvar.find('_Im') > -1:
              updateThese = [0,1,0]
            else:
              print "Error: str=" + str(mutvar)
              updatedThese = [0,0,0]
            for iupdated in range(0,3):
              if updateThese[iupdated]:
                print """forsec L5PC."""+str(updatedVars[iupdated])+""" {
  """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
  }"""
                h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
  """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
  }""")
        print mutText
        tstop = 4000.0
        squareAmp = 1.2
        squareDur = 3800.0
        epsp_Imax = 0.0
        h("""
  tstop = """+str(tstop)+"""
  v_init = """+str(v0)+"""
  cai0_ca_ion = """+str(ca0)+"""
  st1.amp = """+str(squareAmp)+"""
  st1.del = 200
  st1.dur = """+str(squareDur)+"""
  syn1.imax = """+str(epsp_Imax)+"""
  """)
        h.init()
        h.run()
  
        times=np.array(h.tvec)
        Vsoma=np.array(h.vsoma)
        Vdend=np.array(h.vdend)
        Casoma=np.array(h.casoma)
        Cadend=np.array(h.cadend)
        spikes = mytools.spike_times(times,Vsoma,-20,-45)
        spTimesThisCoeff = spikes[:]
        nSpikes1 = len(spikes)
  
        if nSpikes1 > 5:
          spts = spikes[len(spikes)-3:len(spikes)]
          istart = next((i for i,x in enumerate(times) if x > spts[0]))
          iend = next((i for i,x in enumerate(times) if x > spts[1]))+4
          nsteps = iend-istart-1
          tdiff = [y-x for x,y in zip(times[istart:iend-1],times[istart+1:iend])]
          cadiff = [y-x for x,y in zip(Casoma[istart:iend-1],Casoma[istart+1:iend])]
          caddiff = [y-x for x,y in zip(Cadend[istart:iend-1],Cadend[istart+1:iend])]
          caderiv1 = [y/x for x,y in zip(tdiff[0:nsteps-1],cadiff[0:nsteps-1])]
          caderiv2 = [y/x for x,y in zip(tdiff[1:nsteps],cadiff[1:nsteps])]
          caderiv = [(x+y)/2.0 for x,y in zip(caderiv1,caderiv2)]
          cadderiv = [y/x for x,y in zip(tdiff,caddiff)]
  
        #Print the parameters and their default values:
        for idefval in range(0,len(defVals.keys())):
          thisdefval = defVals.keys()[idefval]
          if thisdefval.find('_Im') > -1:
            h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][1]))
            #) #+" (def="+str(defVals[thisdefval])+")"
          else:
            h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][0]))
            #h('print L5PC.soma[0]."+thisdefval) #+" (def="+str(defVals[thisdefval])+")"                       
  
        #Restore default values:
        for imutvar in range(0,len(MT[igene][imut])):
          mutvars = allmutvars[iallmutval][imutvar]
          mutvals = allmutvals[iallmutval][imutvar]
          if type(mutvars) is str:
            mutvars = [mutvars]
          for kmutvar in range(0,len(mutvars)):
            mutvar = mutvars[kmutvar]
            newVal = defVals[mutvar]
            if mutvar.find('_Ih') > -1:
              updateThese = [1,1,1]
            elif mutvar.find('_Ca_HVA') > -1 or mutvar.find('_Ca_LVAst') > -1 or mutvar.find('_SKv3.1') > -1 or mutvar.find('_Ca_HVA') > -1 or mutvar.find('_SK_E2') > -1 or mutvar.find('_NaTa_t') > -1 or mutvar.find('_CaDynamics_E2') > -1:
              updateThese = [1,1,0]
            elif mutvar.find('_K_Pst') > -1 or mutvar.find('_K_Tst') > -1 or mutvar.find('_Nap_Et2') > -1:
              updateThese = [1,0,0]
            elif mutvar.find('_Im') > -1:
              updateThese = [0,1,0]
            else:
              print "Error: str=" + str(mutvar)
              updatedThese = [0,0,0]
            for iupdated in range(0,3):
              if updateThese[iupdated]:
                print """forsec L5PC."""+str(updatedVars[iupdated])+""" {
  """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
  }"""
                h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
  """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
  }""")
        spTimesThisMutVal.append(spTimesThisCoeff[:])
        timesThisMutVal.append(times[:])
        VsomasThisMutVal.append(Vsoma[:])
        CasomasThisMutVal.append(Casoma[:])
        VdendsThisMutVal.append(Vdend[:])
        CadendsThisMutVal.append(Cadend[:])
        if iter==-1:
          picklelist = [theseCoeffsAllAll,times,Vsoma,Casoma,Vdend,Cadend,MT]
          file = open('steadystate_cs'+str(icell)+'_control.sav', 'w')
          pickle.dump(picklelist,file)
          file.close()
  
      spTimesThisMut.append(spTimesThisMutVal[:])
      timesThisMut.append(timesThisMutVal[:])
      VsomasThisMut.append(VsomasThisMutVal[:])
      CasomasThisMut.append(CasomasThisMutVal[:])
      VdendsThisMut.append(VdendsThisMutVal[:])
      CadendsThisMut.append(CadendsThisMutVal[:])
      picklelist = [theseCoeffsAllAll,timesThisMutVal,VsomasThisMutVal,CasomasThisMutVal,VdendsThisMutVal,CadendsThisMutVal,MT]
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
  
  
  spTimesAllAll.append(spTimesAll[:])
  timesAllAll.append(timesAll[:])
  VsomasAllAll.append(VsomasAll[:])
  CasomasAllAll.append(CasomasAll[:])
  VdendsAllAll.append(VdendsAll[:])
  CadendsAllAll.append(CadendsAll[:])
#picklelist = [theseCoeffsAllAll,VsomasAllAll,CasomasAllAll,VdendsAllAll,CadendsAllAll,MT]
#file = open('steadystate.sav', 'w')
#pickle.dump(picklelist,file)
#file.close()
  

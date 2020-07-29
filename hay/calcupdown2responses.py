#Changed from _old4: Made high resolution (smaller CVode tolerance) 26.10.2016. Only do cells 0 and 1.
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
#proximalpoints = [100,100,100,200,200,200,300,300,300,400,400,400]
#distalpoints =   [600,750,900,600,750,900,600,750,900,600,750,900]
proximalpoints = [200,200]
distalpoints =   [600,850]
BACdt = 5.0
fs = 8
tstop = 5000.0
#epspdts = [0.25*x for x in range(-80,81)]
#epspdts = [0.5*x for x in range(-40,41)]
epspdts = [4.0*x for x in range(-40,-20)]+[2.0*x for x in range(-40,-20)]+[1.0*x for x in range(-40,-20)]+[0.5*x for x in range(-40,41)]+[1.0*x for x in range(21,41)]+[2.0*x for x in range(21,41)]+[4.0*x for x in range(21,41)]
epspdts_savetimecourses = [-40,-20,-15,-10,-5,0,5,10,15,20,40]

Is_st2 = 1.32
st2coeff = 0.40      #Somatic 5ms pulse
st2coeff_down = 1.35 #Somatic 5ms pulse
st1coeff = 0.85      #Proximal apical 200ms pulse

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
for idefval in range(0,len(keyList)):
  if type(defVals[keyList[idefval]]) is not list:
    defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments
unpicklefile = open('scalings_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
theseCoeffsAllAll = unpickledlist[0]
theseMutValsAllAll = unpickledlist[2]

paramdicts = []
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 1.0, 'S_gCa_HVAbar_Ca_HVA': 1.0})   # 1 spike per burst, control
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 1.6})                               # 1-2 spikes per burst
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2})                               # 2-3 spikes per burst
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.9})   # 3-4 spikes per burst
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.625}) # 3-5 spikes per burst
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.5})   # 4-6 spikes per burst
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.3})   # 5-9 spikes per burst

VsomaupAllAll = []
VsomadownAllAll = []
VdendupAllAll = []
VdenddownAllAll = []

counter = -1
for icell in range(0,7):
  syn1coeff = 0.5 - 0.35*(icell > 1)    #Synaptic epsp-like input 
  morphology_file = "morphologies/cell1.asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  theseCoeffsAll = theseCoeffsAllAll[icell]

  for idist in range(0,len(proximalpoints)):
    counter = counter + 1
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
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.00001)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
access L5PC.soma
objref st1, st2
st2 = new IClamp(0.5)
L5PC.soma st2
st2.amp = """+str(st2coeff*Is_st2)+"""
st2.del = 1000
st2.dur = 5
objref vsoma, casoma, sksoma, vdend, vdend2, vdend3, cadend3, skdend3, tvec, isyn
vsoma = new Vector()
casoma = new Vector()
sksoma = new Vector()
vdend = new Vector()
vdend2 = new Vector()
vdend3 = new Vector()
cadend3 = new Vector()
skdend3 = new Vector()
objref sl,ns,syn1,con1
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
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)
L5PC.apic[siteVec[0]] {
  syn1 = new epsp(siteVec[1])
  syn1.tau0 = 0.5
  syn1.tau1 = 5
  syn1.onset = 1000 + """+str(BACdt)+"""
  syn1.imax = """+str(syn1coeff*Is_syn1)+"""
  cvode.record(&syn1.i,isyn,tvec)
}
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
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend2,tvec)
L5PC.apic[siteVec[0]] st1 = new IClamp(siteVec[1])
st1.amp = """+str(st1coeff*Is_st1)+"""
st1.del = 700
st1.dur = 600

sl = L5PC.locateSites("apic","""+str(fixedpoint)+""")
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
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend3,tvec)
L5PC.apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend3,tvec)
L5PC.apic[siteVec[0]] cvode.record(&ik_SK_E2(siteVec[1]),skdend3,tvec)

L5PC.soma cvode.record(&v(0.5),vsoma,tvec)
L5PC.soma cvode.record(&cai(0.5),casoma,tvec)
L5PC.soma cvode.record(&ik_SK_E2(0.5),sksoma,tvec)
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
        if exists('updown2response_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav'):
          print 'updown2response_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav exists'
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
          thisCa = h.L5PC.soma[0].minCai_CaDynamics_E2

          for idt in range(0,len(epspdts)):
            for iup in range(0,2):
              syn1coeff = 0.5 - 0.35*(icell > 1 or (icell > 0 and iup == 1))    #Synaptic epsp-like input 

              print "st1.amp = "+str(st1coeff*Is_st1*(1-iup))
              print "st2.amp = "+str((st2coeff*(1-iup) + st2coeff_down*iup)*Is_st2)
              print "syn1.imax = "+str(syn1coeff*Is_syn1)
              h("st1.amp = "+str(st1coeff*Is_st1*(1-iup)))
              h("st2.amp = "+str((st2coeff*(1-iup) + st2coeff_down*iup)*Is_st2))
              h("syn1.imax = "+str(syn1coeff*Is_syn1))
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
                file = open('updown2responsetimecourse_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'w')
                pickle.dump(picklelist,file)
                file.close()

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
          file = open('updown2responsetimecourse_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'w')
          pickle.dump(picklelist,file)
          file.close()

        picklelist = [theseCoeffsAllAll,CasomaupThisMutVal,SKsomaupThisMutVal,VdendupThisMutVal,Vdend2upThisMutVal,Vdend3upThisMutVal,Cadend3upThisMutVal,SKdend3upThisMutVal,
                      CasomadownThisMutVal,SKsomadownThisMutVal,VdenddownThisMutVal,Vdend2downThisMutVal,Vdend3downThisMutVal,Cadend3downThisMutVal,SKdend3downThisMutVal,
                      spikesupThisMutVal,spikesdownThisMutVal,epspdts,MT]
        file = open('updown2response_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'w')
        pickle.dump(picklelist,file)
        file.close()

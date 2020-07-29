#coding_nonprop: stimuli not proportional to the threshold current
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
import itertools

random.seed(1)

v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
fs = 8
ITERS = 20
tstop = 11000.0

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


#lensToStart = [100.0 + x*50 for x in range(0,16)]
lensToStart = [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0]
lensToEnd = [200.0, 400.0, 600.0, 800.0, 1000.0, 1300.0]
lenToStartBasal = 0.0
lenToEndBasal = 282.0


gCoeffsAllAllAll = []
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
  for istartdist in range(0,len(lensToStart)):
    lenToStart = lensToStart[istartdist]
    lenToEnd = lensToEnd[istartdist]
    unpicklefile = open('synlocs'+str(lenToStart)+'-'+str(lenToEnd)+'.sav', 'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    Nsyns = unpickledlist[0]
    synlocsAll = unpickledlist[3]

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
  morphology_file = "morphologies/cell1.asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  theseCoeffsAll = theseCoeffsAllAll[icell]
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
objref vsoma, vdend, recSite, vdends, isoma, cadends, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadends = new List()
vdends = new List()
objref sl,st2,ns,syn1,con1,isyn, tvec, syns["""+str(7*Nsyns)+"""]
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
access L5PC.apic[siteVec[0]]
st2 = new IClamp(siteVec[1])
st2.amp = 0
L5PC.apic[siteVec[0]] {
  st2
  syn1 = new epsp(siteVec[1])
  syn1.tau0 = 0.5
  syn1.imax = 0
  syn1.tau1 = 5
  syn1.onset = 145 
  cvode.record(&syn1.i,isyn,tvec)
}
access L5PC.soma
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
access L5PC.apic[siteVec[0]]
cvode.record(&v(siteVec[1]),vdend,tvec)
""")
  for idist in range(0,6):
    h("""
cadends.append(new Vector())
vdends.append(new Vector())
sl = L5PC.locateSites("apic","""+str((lensToStart[idist]+lensToEnd[idist])*0.5)+""")
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
cvode.record(&cai(siteVec[1]),cadends.o[cadends.count()-1],tvec)
cvode.record(&cai(siteVec[1]),vdends.o[vdends.count()-1],tvec)
""")
    for istim in range(0,Nsyns):
      h("""
siteVec[0] = """+str(synlocsAllAll[idist][0][istim][0])+"""
siteVec[1] = """+str(synlocsAllAll[idist][0][istim][1])+"""
access L5PC.apic[siteVec[0]]
L5PC.apic[siteVec[0]] {
  syns["""+str(Nsyns*idist+istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(Nsyns*idist+istim)+"""].e = 0
  syns["""+str(Nsyns*idist+istim)+"""].tau = 5
  syns["""+str(Nsyns*idist+istim)+"""].onset = 10000
  syns["""+str(Nsyns*idist+istim)+"""].gmax = """+str(synconductance)+"""
}
""")
  for istim in range(0,Nsyns):
    h("""
siteVec[0] = """+str(synlocsAllAll[6][0][istim][0])+"""
siteVec[1] = """+str(synlocsAllAll[6][0][istim][1])+"""
access L5PC.dend[siteVec[0]]
L5PC.dend[siteVec[0]] {
  syns["""+str(Nsyns*6+istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(Nsyns*6+istim)+"""].e = 0
  syns["""+str(Nsyns*6+istim)+"""].tau = 5
  syns["""+str(Nsyns*6+istim)+"""].onset = 10000
  syns["""+str(Nsyns*6+istim)+"""].gmax = """+str(synconductance)+"""
}
""")

  paramdict = paramdicts[icell]
  setparams(paramdict)

  styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
  #cols = ['#aaffaa','#aaffaa','#66ff66','#66aaaa','#00aaaa','#00aaaa']
  #cols = ['#00aaaa','#00bb77','#11cc44','#11dd11','#55ee00','#99dd00']
  cols = ['#00aaaa','#11cc44','#55ee00','#bbaa00','#ee6600','#ff0000', '#aa00aa','#772277','#333333']
  #yplus = [1, 2, 3, 4, 5, 6, -1, -2, -3]
  #yplus = [x+3 for x in yplus]
  yplus = [0, 0, 0, 0, 0, 0, 0, 0, 0]
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
      #for iter in [2, 5, 8, -1]:
      iters = [0, 2, 6, 8, -1]
      nspsThisVal = []
      mytime = time.time()
      coding_outputs = []
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
          continue
        #thisg = gcoeffsAllAll[icell][igene][imut][iallmutval][iiter]
        #print "iter="+str(iter)+", thisCoeff="+str(thisCoeff)+", thisg="+str(thisg)
        
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

        h("""
tstop = """+str(tstop)+"""
cai0_ca_ion = """+str(thisCa)+"""
v_init = """+str(v0)+"""
st1.amp = """+str(somaticI)+"""
st1.del = 9900
st1.dur = 200
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
          spikes = mytools.spike_times(times,Vsoma,-20,-45)
          nSpikes1 = len(spikes)
          coding_outputs_thisiter.append([max(array(h.cadends)[i]) for i in range(0,6)]+[nSpikes1]+[max(array(h.vdends)[i]) for i in range(0,6)])
        coding_outputs.append(coding_outputs_thisiter[:])

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

        picklelist = [[],coding_inputs,coding_outputs,[]]
        file = open('codings_nonprop'+str(synconductance)+'_cs'+str(icell)+'_'+str(counter)+'_somaticI'+str(somaticI)+'.sav', 'w')
        pickle.dump(picklelist,file)
        file.close()
      print "codings done in "+str(time.time()-mytime)+" seconds, icell="+str(icell)+", igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)

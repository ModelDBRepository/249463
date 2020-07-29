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
import random
import time


v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
Is = unique([0.34+0.0025*x for x in range(0,11)]+[0.35+0.05*x for x in range(0,22)])
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

import mutation_stuff
MT = mutation_stuff.getMT()
updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments
unpicklefile = open('scalings_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()

theseCoeffsAllAll = unpickledlist[0]
theseMutValsAll = unpickledlist[2]
ITERS = 20
tstop = 11000.0
ISIs = unique([4*x for x in range(0,25)]+[20*x for x in range(0,26)])
currCoeff = 1.1

paramdicts = []
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 1.0, 'S_gCa_HVAbar_Ca_HVA': 1.0})   # 1 spike per burst, control
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 1.6})                               # 1-2 spikes per burst
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2})                               # 2-3 spikes per burst
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.9})   # 3-4 spikes per burst
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.625}) # 3-5 spikes per burst
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.5})   # 4-6 spikes per burst
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.3})   # 5-9 spikes per burst

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

#lensToStart = [100.0 + x*50 for x in range(0,16)]
lensToStart = [150.0, 300.0, 450.0, 600.0, 650.0]

gCoeffsAllAllAll = []


for istartdist in range(0,len(lensToStart)):
 startdist = lensToStart[istartdist]
 gCoeffsAllAll = []
 if len(sys.argv) > 2 and int(float(sys.argv[2])) != istartdist:
   continue
 unpicklefile = open('synlocs'+str(startdist)+'.sav', 'r')
 unpickledlist = pickle.load(unpicklefile)
 unpicklefile.close()
 Nsyns = unpickledlist[0]
 synlocsAll = unpickledlist[3]
 startdist = int(startdist)

 maxLens = [1300,1185]

 for icell in range(0,7):
  synlocs = synlocsAll[0]
  gCoeffsThisVal = []
  if len(sys.argv) > 3 and int(float(sys.argv[3])) != icell:
    continue
  morphology_file = "morphologies/cell1.asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  unpicklefile = open('thresholddistalamp'+str(startdist)+'_cs'+str(icell)+'_comb'+str(icomb)+'.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  gs_thiscomb = unpickledlist[1]

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
forsec L5PC.somatic {
}
forsec L5PC.apical {
}
L5PC.distribute_channels("apic","gIhbar_Ih",2,-0.8696,3.6161,0.0,1.0*2.0870,0.0002)
L5PC.distribute_channels("apic","gCa_HVAbar_Ca_HVA",3,1.0,0.1,685.0,885.0,1.0*0.000555)
L5PC.distribute_channels("apic","gCa_LVAstbar_Ca_LVAst",3,1.0,0.01,685.0,885.0,1.0*0.0187)
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
objref sl,st2,ns,syn1,con1,isyn, tvec, syns["""+str(2*Nsyns)+"""]
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
cvode.record(&cai(siteVec[1]),cadend,tvec)
""")
  for istim in range(0,Nsyns):
    h("""
siteVec[0] = """+str(synlocs[istim][0])+"""
siteVec[1] = """+str(synlocs[istim][1])+"""
access L5PC.apic[siteVec[0]]
L5PC.apic[siteVec[0]] {
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

  paramdict = paramdicts[icell]
  setparams(paramdict)

  styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
  #cols = ['#00aaaa','#11cc44','#55ee00','#bbaa00','#ee6600','#ff0000', '#aa00aa','#772277','#333333']
  cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
  
  defVals = mutation_stuff.getdefvals()
  keyList = defVals.keys()
  for idefval in range(0,len(keyList)):
    if type(defVals[keyList[idefval]]) is not list:
      defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]

  iters = [0, 2, 5, 6, 8, -1]
  for iiter in range(0,len(iters)):
    iter = iters[iiter]
    gCoeffsThisIter = []
    thisg = gs_thiscomb[iiter]
    if iter == 5:
      continue
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
          nextCoeffs = [0.0,2.0,1.0]
  
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
              if mutvar.find('offm') > -1 or mutvar.find('offh') > -1 or mutvar.find('ehcn') > -1:
                newVal =  [x+mutvals*thisCoeff for x in defVals[mutvar]]
                defVals[mutvar] = [x+mutvals*thisCoeff for x in defVals[mutvar]]
                if mutvals >= 0 and kmutvar==0:
                  mutText = mutText + "+" + str(mutvals) +" mV"
                elif kmutvar==0:
                  mutText = mutText  + str(mutvals) +" mV"
              else:
                newVal =  [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
                defVals[mutvar] = [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
                if kmutvar==0:
                  mutText = mutText + "*" + str(mutvals)
              if kmutvar < len(mutvars)-1:
                mutText = mutText + ", "
              if mutvar.find('_Ih') > -1:
                updateThese = [1,1,1]
              elif mutvar.find('_Ca_HVA') > -1 or mutvar.find('_Ca_LVAst') > -1 or mutvar.find('_SKv3.1') > -1 or mutvar.find('_Ca_HVA') > -1 or mutvar.find('_SK_E2') > -1 or mutvar.find\
    ('_NaTa_t') > -1 or mutvar.find('_CaDynamics_E2') > -1:
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
    print "Working on iter = "+str(iter)
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
  st1.del = 200
  st1.dur = 10
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
              spikes = mytools.spike_times(times,Vsoma,-20,-45)
              nSpikes1 = len(spikes)
              print "nextCoeffs="+str(nextCoeffs)+", "+str(nSpikes1)+" spikes, simulation done in "+str(time.time()-timenow)+" seconds"
              nSpikes_normal = 1
              if icell > 0 or startdist <= 200: # For icell=0, 1 spike normally generated (except for the inputs nearest to soma), while for icell=1,2,3,4,5,6, two spikes normally generated
                nSpikes_normal = 2
              hasSpiked = hasSpiked or (nSpikes1 > nSpikes_normal)
              if iterI == 0 and hasSpiked:
                print "istartdist="+str(istartdist)+", icell="+str(icell)+", igene="+str(igene)+", imut="+str(imut)+", iallmuval="+str(iallmutval)+", iiter="+str(iiter)+", iISI="+str(iISI)+": extra spikes elicited for iterI=0!"
              if iterI > 0 and not hasSpiked:
                print "istartdist="+str(istartdist)+", icell="+str(icell)+", igene="+str(igene)+", imut="+str(imut)+", iallmuval="+str(iallmutval)+", iiter="+str(iiter)+", iISI="+str(iISI)+": no extra spikes for iterI>0!"
                nextCoeffs = [nextCoeffs[1],2*nextCoeffs[1],1.5*nextCoeffs[1]]
                continue
              if iterI > 1 and nSpikes1 > nSpikes_normal:
                nextCoeffs = [nextCoeffs[0],nextCoeffs[2],0.5*(nextCoeffs[0]+nextCoeffs[2])]
              if iterI > 1 and nSpikes1 <= nSpikes_normal:
                nextCoeffs = [nextCoeffs[2],nextCoeffs[1],0.5*(nextCoeffs[2]+nextCoeffs[1])]
            gCoeffsThisISI = nextCoeffs[:]
            gCoeffsThisIter.append(gCoeffsThisISI[:])
  
    defVals = mutation_stuff.getdefvals()
    keyList = defVals.keys()
    for idefval in range(0,len(keyList)):
      if type(defVals[keyList[idefval]]) is not list:
        defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
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
    gCoeffsThisVal.append(gCoeffsThisIter[:])

  picklelist = [theseCoeffsAll,gCoeffsThisVal,MT]
  file = open('PPIcoeffs'+str(startdist)+'_cs'+str(icell)+'_comb'+str(icomb)+'.sav', 'w')
  pickle.dump(picklelist,file)
  file.close()

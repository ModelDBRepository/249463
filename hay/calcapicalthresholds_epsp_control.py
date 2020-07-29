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


random.seed(1)

v0 = -80
ca0 = 0.0001
distalpoint = 760
dists = [50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000]
BACdt = 5.0
fs = 8
tstop = 5000.0
#epspdts = range(-20,21)
##epspdts = [-20,0,20]
epspdts = [0.25*x for x in range(-80,81)]

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

IsAllAll = []
counter = -1
for icell in range(0,7):
  morphology_file = "morphologies/cell1.asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  theseCoeffsAll = theseCoeffsAllAll[icell]
  times_control = [[[],[],[]],[[],[],[]]]
  Vsoma_control = [[[],[],[]],[[],[],[]]]
  Vdend_control = [[[],[],[]],[[],[],[]]]

  IsAll = []
  for idist in range(0,len(dists)):
    counter = counter+1
    if len(sys.argv) > 1 and counter < int(sys.argv[1]):
      continue
    h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.0002)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
objref st1, st2
st2 = new IClamp(0.5)
L5PC.soma st2
st2.amp = 0
st2.del = 1000
st2.dur = 5
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
vdend2 = new Vector()
cadend2 = new Vector()
objref sl,ns,syn1,syn2,con1,isyn, tvec
isyn = new Vector()
tvec = new Vector()
sl = new List()
double siteVec[2]
sl = L5PC.locateSites("apic","""+str(dists[idist])+""")
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
cvode.record(&v(siteVec[1]),vdend,tvec)
cvode.record(&cai(siteVec[1]),cadend,tvec)
L5PC.apic[siteVec[0]] {
  syn1 = new epsp(siteVec[1])
  syn1.tau0 = 0.5
  syn1.tau1 = 5
  syn1.onset = 1000
  syn1.imax = 0
  syn2 = new epsp(siteVec[1])
  syn2.tau0 = 0.5
  syn2.tau1 = 5
  syn2.onset = 145 + """+str(BACdt)+""" 
  syn2.imax = 0
  cvode.record(&syn1.i,isyn,tvec)
}
sl = L5PC.locateSites("apic","""+str(dists[idist])+""")
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
L5PC.apic[siteVec[0]] st1 = new IClamp(siteVec[1])
st1.amp = 0
st1.del = 900
st1.dur = 200
access L5PC.soma
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
""")

    paramdict = paramdicts[icell]
    setparams(paramdict)

    Is = []
    maxstim2 = 20.0
    ITERS = 33
    nextIs = [0.0,20.0,10.0]
    for iter in range(0,ITERS):
    
      myI = nextIs[min(iter,2)]
      h("st1.amp = 0")
      h("st2.amp = 0")
      h("syn1.imax = "+str(myI))
      h("""
tstop = """+str(tstop)+"""
cai0_ca_ion = """+str(ca0)+"""
v_init = """+str(v0)+"""
syn1.onset = """+str(1000)+"""
""")
      h.init()
      try:
        h.run()
      except RuntimeError:
        print "Too large I!"
        nextIs = [nextIs[0],nextIs[2],0.5*(nextIs[0]+nextIs[2])]
        continue
      times=np.array(h.tvec)
      Vsoma=np.array(h.vsoma)
      Vdend=np.array(h.vdend)

      nSpikes = len(mytools.spike_times(times,Vsoma,-20,-45))
      if iter > 2 and nSpikes > 0:
        nextIs = [nextIs[0],nextIs[2],0.5*(nextIs[0]+nextIs[2])]
      if iter > 2 and nSpikes == 0:
        nextIs = [nextIs[2],nextIs[1],0.5*(nextIs[2]+nextIs[1])]
      print str(nSpikes)+", nextIs="+str(nextIs)
      Is.append(nextIs[2])

      picklelist = [Is,dists,idist]
      file = open('apicalthresholds_epsp_icell'+str(icell)+'_dist'+str(dists[idist])+'.sav', 'w')
      pickle.dump(picklelist,file)
      file.close()

    IsAll.append(Is[:])
    picklelist = [IsAll,dists]
    file = open('apicalthresholds_epsp_icell'+str(icell)+'.sav', 'w')
    pickle.dump(picklelist,file)
    file.close()

  IsAllAll.append(IsAll[:])
  picklelist = [IsAllAll,dists]
  file = open('apicalthresholds_epsp.sav', 'w')
  pickle.dump(picklelist,file)
  file.close()

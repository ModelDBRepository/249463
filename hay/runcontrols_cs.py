from neuron import h
import matplotlib
import numpy
from pylab import *
import pickle
import protocol

import mutation_stuff
import mytools
defValsMut = mutation_stuff.getdefvals()
defVals = protocol.get_defvals()

spikfreqsAll = []
timescAll = []
VsomacAll = []
VDerivcAll = []
VDcoeffAll = []
VdendcAll = []
VdDcoeffAll = []
VdDerivcAll = []
CasomacAll = []
CaDerivcAll = []
CaDcoeffAll = []
CadendcAll = []
CadDerivcAll = []
CadDcoeffAll = []
times_controlAll = []
Vsoma_controlAll = []
Vdend_controlAll = []
Casoma_controlAll = []
Cadend_controlAll = []

def setparams(params):
  params_copy = params.copy()
  keys = params.keys()
  for ikey in range(0,len(keys)):
    params_copy[keys[ikey]] = params[keys[ikey]]*defVals[keys[ikey]]

  #Apply the new parameter values: 
  IhChanged = False
  CaHVAChanged = False
  CaLVAChanged = False
  aIh = defVals['D_aIh']
  bIh = defVals['D_bIh']
  print "setparams: params_copy = "+str(params_copy)
  for ikey in range(0,len(keys)):
    key = keys[ikey]
    section = []
    if key[0:2] == "A_" or key[0:2] == "B_" or key[0:2] == "S_" or key[0:2] == "*_":
      if key[0:2] == "A_":
        section = "apical"
      if key[0:2] == "B_":
        section = "basal"
      if key[0:2] == "S_":
        section = "somatic"
      if key[0:2] == "*_":
        section = "all"
      if len(section) > 0:
        print("forsec L5PC."+section+" "+key[2:]+" = "+str(params_copy[key]))
        h("forsec L5PC."+section+" "+key[2:]+" = "+str(params_copy[key]))
    elif key[0:2] == "D_":
      if key[2:] == "aIh":
        IhChanged = True
        aIh = params_copy[key]
      if key[2:] == "bIh":
        IhChanged = True
        bIh = params_copy[key]
      if key == "aCa_HVA":
        CaHVAChanged = True
      if key[2:] == "aCa_LVAst":
        CaLVAChanged = True
  if IhChanged:
    print("""L5PC.distribute_channels(\"apic\",\"gIhbar_Ih\",2,"""+str(aIh)+""",3.6161,0.0,"""+str(bIh)+""",0.00020000000)""")
    h("""L5PC.distribute_channels(\"apic\",\"gIhbar_Ih\",2,"""+str(aIh)+""",3.6161,0.0,"""+str(bIh)+""",0.00020000000)""")
  if CaLVAChanged:
    print("""L5PC.distribute_channels(\"apic\",\"gCa_LVAstbar_Ca_LVAst\",3,1.000000,0.01,685.000000,885.000000,"""+str(params_copy['D_aCa_LVAst'])+""")""")
    h("""L5PC.distribute_channels(\"apic\",\"gCa_LVAstbar_Ca_LVAst\",3,1.000000,0.01,685.000000,885.000000,"""+str(params_copy['D_aCa_LVAst'])+""")""")
  if CaHVAChanged:
    print("""L5PC.distribute_channels(\"apic\",\"gCa_HVAbar_Ca_HVA\",3,1.000000,0.100000,685.000000,885.000000,"""+str(params_copy['D_aCa_HVA'])+""")""")
    h("""L5PC.distribute_channels(\"apic\",\"gCa_HVAbar_Ca_HVA\",3,1.000000,0.100000,685.000000,885.000000,"""+str(params_copy['D_aCa_HVA'])+""")""")


paramdicts = []
paramdicts.append({})                                                         # 1 spike per burst, control
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 1.6})                               # 1-2 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2})                               # 2-3 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.9})   # 3-4 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.625}) # 3-5 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.5})   # 4-6 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.3})   # 5-9 spikes per burst      



for icell in range(0,len(paramdicts)):
  morphology_file = "morphologies/cell1.asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"
  v0 = -80
  ca0 = 0.0001

  proximalpoint = 400
  distalpoint = 620
  BACdt = 5.0

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
objref st1
st1 = new IClamp(0.5)
st1.del = 3400
L5PC.soma st1
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
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend,tvec)
L5PC.apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend,tvec)
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

  h("""
print "eca: ", eca
print "somatic gCa_HVA: ", L5PC.soma.gCa_HVAbar_Ca_HVA
print "distalpoint gCa_HVA: ", L5PC.apic[siteVec[0]].gCa_HVAbar_Ca_HVA
print "distalpoint gCa_LVA: ", L5PC.apic[siteVec[0]].gCa_LVAstbar_Ca_LVAst
print "distalpoint gIm: ", L5PC.apic[siteVec[0]].gImbar_Im
print "distalpoint gIh: ", L5PC.apic[siteVec[0]].gIhbar_Ih
print "distalpoint gSKv3_1: ", L5PC.apic[siteVec[0]].gSKv3_1bar_SKv3_1
print "distalpoint gSK_E2: ", L5PC.apic[siteVec[0]].gSK_E2bar_SK_E2
print "distalpoint gNaTa: ", L5PC.apic[siteVec[0]].gNaTa_tbar_NaTa_t
print "proximalpoint gCa_HVA: ", L5PC.apic[siteVec[0]].gCa_HVAbar_Ca_HVA
print "proximalpoint gCa_LVA: ", L5PC.apic[siteVec[0]].gCa_LVAstbar_Ca_LVAst
""")
  #Print the parameters and their default values:
  for idefval in range(0,len(defValsMut.keys())):
    thisdefval = defValsMut.keys()[idefval]
    if thisdefval.find('_Im') > -1:
      if type(defValsMut[thisdefval]) is float:
        h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defValsMut[thisdefval]))
      else:
        h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defValsMut[thisdefval][1]))
    else:
      if type(defValsMut[thisdefval]) is float:
        h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defValsMut[thisdefval]))
      else:
        h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defValsMut[thisdefval][0]))

  #Is = [0.2,0.4,0.6,0.8,1.0,1.2,1.4]
  Is = [0.1*x for x in range(0,16)]
  spikfreqs = len(Is)*[0]
  for iI in range(0,len(Is)):
    squareAmp = Is[iI]
    squareDur = 3800
    tstop = 7200
    h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.dur = """+str(squareDur)+"""
""")
    h.init()
    h.run()

    times=np.array(h.tvec)
    Vsoma=np.array(h.vsoma)
    Vdend=np.array(h.vdend)
    Casoma=np.array(h.casoma)
    Cadend=np.array(h.cadend)
    spikes = mytools.spike_times(times,Vsoma,0,100)
    spikfreqs[iI] = sum([1 for x in spikes if x >= 500.0])/3.5

    f, axarr = subplots(1, 1)
    axarr.plot(times, Vsoma)
    axarr.plot(spikes, [1 for x in spikes],'r.')
    axarr.set_title("Perisomatic firing, nspikes="+str(len(spikes))+", after 300ms nspikes="+str(sum([1 for x in spikes if x >= 3700.0])))
    axarr.set_ylim([-100,40])
    axarr.set_xlim([3390,7200])
    f.savefig("runcontrols_cs"+str(icell)+"_"+str(iI)+".eps")

    if abs(Is[iI]-1.0) < 0.0001:
      Vsoma_control = Vsoma
      Casoma_control = Casoma
      Vdend_control = Vdend
      Cadend_control = Cadend
      times_control = times
      spikes_control = spikes
  

  isis = [y-x for x,y in zip(spikes_control[0:-1],spikes_control[1:])]
  if len(isis) > 0 and mean(isis) != 0:
    CVISI = std(isis)/mean(isis) #std is defined in pylab                                                                                                          
  else:
    CVISI = nan #NaN defined in pylab  
  ispikelastburst = -1
  nspikesperbursts = []
  for ispike in range(0,len(isis)):
    if isis[ispike] > 1.1*mean(isis):
      nspikesperbursts.append(ispike-ispikelastburst)
      ispikelastburst = ispike
  if len(nspikesperbursts) < 2:
    nSpikesPerLc = 1
  else:
    nSpikesPerLc = mean(nspikesperbursts[1:])

  print "nSpikesPerLc = "+str(nSpikesPerLc)
  
  if nSpikesPerLc-int(nSpikesPerLc) > 0:
    nSpikesPerLc = int(nSpikesPerLc)+1
  else:
    nSpikesPerLc = int(nSpikesPerLc)

  #First calculate the VDcoeff from a total limit cycle spanning more than 1 full cycle
  spts = spikes_control[len(spikes_control)-nSpikesPerLc-1:len(spikes_control)]
  istart = next((i for i,x in enumerate(times_control) if x > spts[0]))
  iend = next((i for i,x in enumerate(times_control) if x > spts[len(spts)-1]))+5
  Vsomac = Vsoma_control[istart:iend]
  timesc = times_control[istart:iend]
  VDerivc = mytools.membpotderivs(timesc,Vsomac)
  VDcoeff =  mytools.limitcyclescaledv(Vsomac,VDerivc,Vsomac,VDerivc)

  spdts = [y-x for x,y in zip(spts[0:len(spts)-1],spts[1:len(spts)])]
  imaxspdt = spdts.index(max(spdts))
  iminspdt = spdts.index(min(spdts))
  tstart_lc1 = spts[imaxspdt]
  tstop_lc1 = spts[imaxspdt+1]
  tstart_lc2 = spts[iminspdt]
  tstop_lc2 = spts[iminspdt+1]

  #Then calculate the limit cycle corresponding to the longest ISI
  istart = next((i for i,x in enumerate(times_control) if x > tstart_lc1))
  iend = next((i for i,x in enumerate(times_control) if x > tstop_lc1))+5
  nsteps = iend-istart-1
  Vsomac1 = Vsoma_control[istart:iend]
  Casomac1 = Casoma_control[istart:iend]
  Vdendc1 = Vdend_control[istart:iend]
  Cadendc1 = Cadend_control[istart:iend]
  timesc1 = times_control[istart:iend]
  VDerivc1 = mytools.membpotderivs(timesc1,Vsomac1)
  CaDerivc1 = mytools.membpotderivs(timesc1,Casomac1)
  VdDerivc1 = mytools.membpotderivs(timesc1,Vdendc1)
  CadDerivc1 = mytools.membpotderivs(timesc1,Cadendc1)
  VDcoeff1 =  mytools.limitcyclescaledv(Vsomac1,VDerivc1,Vsomac1,VDerivc1)
  CaDcoeff1 =  mytools.limitcyclescaledv(Casomac1,CaDerivc1,Casomac1,CaDerivc1)
  VdDcoeff1 =  mytools.limitcyclescaledv(Vdendc1,VdDerivc1,Vdendc1,VdDerivc1)
  CadDcoeff1 =  mytools.limitcyclescaledv(Cadendc1,CadDerivc1,Cadendc1,CadDerivc1)

  close("all")
  f, axarr = subplots(3, 2)
  axarr[0,0].plot(Vsomac1[1:nsteps-1],VDerivc1)
  axarr[0,0].set_title("Vsoma")
  axarr[0,1].plot(Casomac1[1:nsteps-1],CaDerivc1)
  axarr[0,1].set_title("Casoma")
  axarr[1,0].plot(Vdendc1[1:nsteps-1],VdDerivc1)
  axarr[1,0].set_title("Vdend")
  axarr[1,1].plot(Cadendc1[1:nsteps-1],CadDerivc1)
  axarr[1,1].set_title("Cadend")
  axarr[2,0].plot(times_control,Vsoma_control)
  axarr[2,0].plot(times_control[istart:iend],Vsoma_control[istart:iend])

  Vsomac1 = Vsomac1[1:nsteps-1]
  Vdendc1 = Vdendc1[1:nsteps-1]
  Casomac1 = Casomac1[1:nsteps-1]
  Cadendc1 = Cadendc1[1:nsteps-1]
  timesc1 = timesc1[1:nsteps-1]

  #Then calculate the limit cycle corresponding to the shortest ISI
  istart = next((i for i,x in enumerate(times_control) if x > tstart_lc2))
  iend = next((i for i,x in enumerate(times_control) if x > tstop_lc2))+5
  nsteps = iend-istart-1

  Vsomac2 = Vsoma_control[istart:iend]
  Casomac2 = Casoma_control[istart:iend]
  Vdendc2 = Vdend_control[istart:iend]
  Cadendc2 = Cadend_control[istart:iend]
  timesc2 = times_control[istart:iend]
  VDerivc2 = mytools.membpotderivs(timesc2,Vsomac2)
  CaDerivc2 = mytools.membpotderivs(timesc2,Casomac2)
  VdDerivc2 = mytools.membpotderivs(timesc2,Vdendc2)
  CadDerivc2 = mytools.membpotderivs(timesc2,Cadendc2)
  VDcoeff2 =  mytools.limitcyclescaledv(Vsomac2,VDerivc2,Vsomac2,VDerivc2)
  CaDcoeff2 =  mytools.limitcyclescaledv(Casomac2,CaDerivc2,Casomac2,CaDerivc2)
  VdDcoeff2 =  mytools.limitcyclescaledv(Vdendc2,VdDerivc2,Vdendc2,VdDerivc2)
  CadDcoeff2 =  mytools.limitcyclescaledv(Cadendc2,CadDerivc2,Cadendc2,CadDerivc2)

  axarr[0,0].plot(Vsomac2[1:nsteps-1],VDerivc2)
  axarr[0,0].set_title("Vsoma")
  axarr[0,1].plot(Casomac2[1:nsteps-1],CaDerivc2)
  axarr[0,1].set_title("Casoma")
  axarr[1,0].plot(Vdendc2[1:nsteps-1],VdDerivc2)
  axarr[1,0].set_title("Vdend")
  axarr[1,1].plot(Cadendc2[1:nsteps-1],CadDerivc2)
  axarr[1,1].set_title("Cadend")
  axarr[2,0].plot(times_control[istart:iend],Vsoma_control[istart:iend])
  axarr[2,0].set_xlim([min(tstart_lc1,tstart_lc2)-10,max(tstop_lc1,tstop_lc2)+10])
  axarr[2,0].set_title("Vsoma")
  f.savefig("runcontrols_cs"+str(icell)+"_cycle.eps")

  Vsomac2 = Vsomac2[1:nsteps-1]
  Vdendc2 = Vdendc2[1:nsteps-1]
  Casomac2 = Casomac2[1:nsteps-1]
  Cadendc2 = Cadendc2[1:nsteps-1]
  timesc2 = timesc2[1:nsteps-1]

  VDcoeffAll.append(VDcoeff)

  timescAll.append(timesc1[:])
  VsomacAll.append(Vsomac1[:])
  VDerivcAll.append(VDerivc1[:])
  VDcoeffAll.append(VDcoeff1)
  VdendcAll.append(Vdendc1[:])
  VdDerivcAll.append(VdDerivc1[:])
  VdDcoeffAll.append(VdDcoeff1)
  CasomacAll.append(Casomac1[:])
  CaDerivcAll.append(CaDerivc1[:])
  CaDcoeffAll.append(CaDcoeff1)
  CadendcAll.append(Cadendc1[:])
  CadDerivcAll.append(CadDerivc1[:])
  CadDcoeffAll.append(CadDcoeff1)

  spikfreqsAll.append(spikfreqs[:])
  timescAll.append(timesc2[:])
  VsomacAll.append(Vsomac2[:])
  VDerivcAll.append(VDerivc2[:])
  VDcoeffAll.append(VDcoeff2)
  VdendcAll.append(Vdendc2[:])
  VdDerivcAll.append(VdDerivc2[:])
  VdDcoeffAll.append(VdDcoeff2)
  CasomacAll.append(Casomac2[:])
  CaDerivcAll.append(CaDerivc2[:])
  CaDcoeffAll.append(CaDcoeff2)
  CadendcAll.append(Cadendc2[:])
  CadDerivcAll.append(CadDerivc2[:])
  CadDcoeffAll.append(CadDcoeff2)

  times_controlAll.append(times_control[:])
  Vsoma_controlAll.append(Vsoma_control[:])
  Vdend_controlAll.append(Vdend_control[:])
  Casoma_controlAll.append(Casoma_control[:])
  Cadend_controlAll.append(Cadend_control[:])

  picklelist = [spikfreqs,timescAll,VsomacAll,VDerivcAll,VDcoeffAll,VdendcAll,VdDerivcAll,VdDcoeffAll,CasomacAll,CaDerivcAll,CaDcoeffAll,
                CadendcAll,CadDerivcAll,CadDcoeffAll,times_control,Vsoma_control,Vdend_control,Casoma_control,Cadend_control,Is]
  file = open('control_cs'+str(icell)+'.sav', 'w')
  pickle.dump(picklelist,file)
  file.close()


picklelist = [spikfreqsAll,timescAll,VsomacAll,VDerivcAll,VDcoeffAll,VdendcAll,VdDerivcAll,VdDcoeffAll,CasomacAll,CaDerivcAll,CaDcoeffAll,
              CadendcAll,CadDerivcAll,CadDcoeffAll,times_controlAll,Vsoma_controlAll,Vdend_controlAll,Casoma_controlAll,Cadend_controlAll,Is]
file = open('control_cs.sav', 'w')
pickle.dump(picklelist,file)
file.close()


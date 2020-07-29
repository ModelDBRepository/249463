from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import sys


v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
Is = [0.3+0.1*x for x in range(0,11)]

spTimesAllAll = []
spTimesAllAll2 = []
nSpikesAllAll = []
ISIs_allAll = []

gNaTa_apics = [1.0+x*0.1 for x in range(0,13)] + [2.2 for x in range(0,31)]
gCaHVA_somas = [1.0 for x in range(0,13)] + [0.975-0.025*x for x in range(0,31)]

condSuffixes = ['Ca_HVA','Ca_LVAst','Ih','K_Pst','K_Tst','NaTa_t','Nap_Et2','SK_E2','SKv3_1']
condNames = ['gCa_HVAbar_Ca_HVA', 'gCa_LVAstbar_Ca_LVAst', 'gIhbar_Ih', 'gImbar_Im', 'gK_Pstbar_K_Pst', 'gK_Tstbar_K_Tst', 'gNaTa_tbar_NaTa_t', 'gNap_Et2bar_Nap_Et2', 'gSK_E2bar_SK_E2', 'gSKv3_1bar_SKv3_1']

for icell in range(0,1):
  spTimesAll = []
  spTimesAll2 = []
  nSpikesAll = []
  spikesPerBurstsAll = []

  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
cvode.atol(0.001)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
objref st1
st1 = new IClamp(0.5)
L5PC.soma st1
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
  syn1.onset = 200 + """+str(BACdt)+""" 
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


  styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
  #cols = ['#00aaaa','#11cc44','#55ee00','#bbaa00','#ee6600','#ff0000', '#aa00aa','#772277','#333333']
  cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
  
  counter = -1
  for ig in range(0,len(gNaTa_apics)):
    counter = counter+1
    if len(sys.argv) > 1 and int(float(sys.argv[1])) != counter:
      continue
    
    spTimesThisCoeff = []
    spTimesThisCoeff2 = []
    nSpikesThisCoeff = []
    spikesPerBurstsThisCoeff = []


    print("""forsec L5PC.apical gNaTa_tbar_NaTa_t = gNaTa_tbar_NaTa_t*"""+str(gNaTa_apics[ig]))
    h("""forsec L5PC.apical gNaTa_tbar_NaTa_t = gNaTa_tbar_NaTa_t*"""+str(gNaTa_apics[ig]))
    print("""forsec L5PC.somatic gCa_HVAbar_Ca_HVA = gCa_HVAbar_Ca_HVA*"""+str(gCaHVA_somas[ig]))
    h("""forsec L5PC.somatic gCa_HVAbar_Ca_HVA = gCa_HVAbar_Ca_HVA*"""+str(gCaHVA_somas[ig]))

    for iI in range(0,len(Is)):
      tstop = 4000.0
      squareAmp = Is[iI]
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
  syn1.onset = 200 + """+str(BACdt)+""" 
  """)
      h.init()
      h.run()

      times=np.array(h.tvec)
      Vsoma=np.array(h.vsoma)
      close("all")
      f,axarr = subplots(2,1)
      axarr[0].plot(times,Vsoma)
      axarr[1].plot(times,Vsoma)
      spikes = mytools.spike_times(times,Vsoma,-20,-45)
      spikes2 = mytools.spike_times(times,Vsoma,-35,100)
      isis = [y-x for x,y in zip(spikes[0:-1],spikes[1:])]
      if len(isis) > 0 and mean(isis) != 0:
        CVISI = std(isis)/mean(isis) #std is defined in pylab
      else:
        CVISI = nan #NaN defined in pylab

      ispikelastburst = -1
      nspikesperbursts = []
      for ispike in range(0,len(isis)):
        if isis[ispike] > 40:
          nspikesperbursts.append(ispike-ispikelastburst)
          ispikelastburst = ispike
      if len(nspikesperbursts) >= 2:
        nspikesperbursts = nspikesperbursts[1:]

      nSpikes2 = sum([1 for x in spikes if x >= 500.0])
      spTimesThisCoeff.append(spikes[:])
      spTimesThisCoeff2.append(spikes2[:])
      nSpikesThisCoeff.append(nSpikes2)
      spikesPerBurstsThisCoeff.append(nspikesperbursts[:])
      axarr[0].plot(spikes,[1 for x in spikes],'ro')
      axarr[1].plot(spikes,[1 for x in spikes],'ro')
      axarr[0].set_title("gNaTa_tbar_NaTa_t*"""+str(gNaTa_apics[ig])+", gCa_HVAbar_Ca_HVA*"""+str(gCaHVA_somas[ig]))
      axarr[1].set_xlim([3700,4000])
      f.savefig('nspikesperbursts_cs'+str(icell)+'_'+str(ig)+'_'+str(iI)+'.eps')

    picklelist = [spikesPerBurstsThisCoeff, spTimesThisCoeff, spTimesThisCoeff2, nSpikesThisCoeff]
    file = open('spikesperburst2_'+str(ig)+'.sav', 'w')
    pickle.dump(picklelist,file)
    file.close()

    spTimesAll.append(spTimesThisCoeff[:])
    spTimesAll2.append(spTimesThisCoeff2[:])
    nSpikesAll.append(nSpikesThisCoeff[:])
    spikesPerBurstsAll.append(nspikesperbursts[:])

    print("""forsec L5PC.apical gNaTa_tbar_NaTa_t = gNaTa_tbar_NaTa_t/"""+str(gNaTa_apics[ig]))
    h("""forsec L5PC.apical gNaTa_tbar_NaTa_t = gNaTa_tbar_NaTa_t/"""+str(gNaTa_apics[ig]))
    print("""forsec L5PC.somatic gCa_HVAbar_Ca_HVA = gCa_HVAbar_Ca_HVA/"""+str(gCaHVA_somas[ig]))
    h("""forsec L5PC.somatic gCa_HVAbar_Ca_HVA = gCa_HVAbar_Ca_HVA/"""+str(gCaHVA_somas[ig]))

  spTimesAllAll.append(spTimesAll[:])
  spTimesAllAll2.append(spTimesAll2[:])
  nSpikesAllAll.append(nSpikesAll[:])

#picklelist = [nSpikesAllAll,spTimesAllAll,spTimesAllAll2,condSuffixes,coeffs,Is]
#file = open('ifcurves.sav', 'w')
#pickle.dump(picklelist,file)
#file.close()
  

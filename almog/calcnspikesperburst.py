from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import sys


v0 = -62
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
Is = [0.65+0.025*x for x in range(0,11)]

spTimesAllAll = []
spTimesAllAll2 = []
nSpikesAllAll = []
ISIs_allAll = []

gsk_apics = [1.0+x*0.25 for x in range(0,13)] + [4.0 for x in range(0,20)]
gcas = [1.0 for x in range(0,13)] + [1.05+0.05*x for x in range(0,20)]

condSuffixes = ['bk','sk','cah','car','iH','iA','kslow','na']
gNames = ['gbar','gbar','pbar','pbar','gbar','gbar','gbar','gbar']


for icell in range(0,1):
  spTimesAll = []
  spTimesAll2 = []
  nSpikesAll = []
  spikesPerBurstsAll = []
  spikesPerBurstsAll2 = []

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
tstop = 8000
""")


  counter = -1
  for ig in range(0,len(gsk_apics)):
    counter = counter+1
    if len(sys.argv) > 1 and int(float(sys.argv[1])) != counter:
      continue
    
    spTimesThisCoeff = []
    spTimesThisCoeff2 = []
    nSpikesThisCoeff = []
    spikesPerBurstsThisCoeff = []
    spikesPerBurstsThisCoeff2 = []

    print("""forall if(ismembrane(\"sk\")) gbar_sk = gbar_sk*"""+str(gsk_apics[ig]))
    h("""forall if(ismembrane(\"sk\")) gbar_sk = gbar_sk*"""+str(gsk_apics[ig]))
    print("""forall if(ismembrane(\"cah\")) pbar_cah = pbar_cah*"""+str(gcas[ig]))
    h("""forall if(ismembrane(\"cah\")) pbar_cah = pbar_cah*"""+str(gcas[ig]))
    print("""forall if(ismembrane(\"car\")) pbar_car = pbar_car*"""+str(gcas[ig]))
    h("""forall if(ismembrane(\"car\")) pbar_car = pbar_car*"""+str(gcas[ig]))

    styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
    cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']

    for iI in range(0,len(Is)):
      tstop = 8000.0
      squareAmp = Is[iI]
      squareDur = 7800.0
      h("""
  tstop = """+str(tstop)+"""
  v_init = """+str(v0)+"""
  cai0_ca_ion = """+str(ca0)+"""
  st1.amp = """+str(squareAmp)+"""
  st1.del = 200
  st1.dur = """+str(squareDur)+"""
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
      nspikesperbursts2 = []
      for ispike in range(0,len(isis)):
        if isis[ispike] > 1.1*mean(isis):
          nspikesperbursts.append(ispike-ispikelastburst)
          ispikelastburst = ispike
      if len(nspikesperbursts) >= 2:
        nspikesperbursts = nspikesperbursts[1:]
      for ispike in range(0,len(isis)):
        if isis[ispike] > 40:
          nspikesperbursts2.append(ispike-ispikelastburst)
          ispikelastburst = ispike
      if len(nspikesperbursts2) >= 2:
        nspikesperbursts2 = nspikesperbursts2[1:]


      nSpikes2 = sum([1 for x in spikes if x >= 500.0])
      spTimesThisCoeff.append(spikes[:])
      spTimesThisCoeff2.append(spikes2[:])
      nSpikesThisCoeff.append(nSpikes2)
      spikesPerBurstsThisCoeff.append(nspikesperbursts[:])
      spikesPerBurstsThisCoeff2.append(nspikesperbursts2[:])
      axarr[0].plot(spikes,[1 for x in spikes],'ro')
      axarr[1].plot(spikes,[1 for x in spikes],'ro')
      axarr[0].set_title("gbar_sk*"""+str(gsk_apics[ig])+", [gcah,gcar]*"""+str(gcas[ig]))
      axarr[1].set_xlim([3500,4000])
      f.savefig('nspikesperbursts_cs'+str(icell)+'_'+str(ig)+'_'+str(iI)+'.eps')

    picklelist = [spikesPerBurstsThisCoeff, spikesPerBurstsThisCoeff2, spTimesThisCoeff, spTimesThisCoeff2, nSpikesThisCoeff]
    file = open('spikesperburst_'+str(ig)+'.sav', 'w')
    pickle.dump(picklelist,file)
    file.close()

    spTimesAll.append(spTimesThisCoeff[:])
    spTimesAll2.append(spTimesThisCoeff2[:])
    nSpikesAll.append(nSpikesThisCoeff[:])
    spikesPerBurstsAll.append(nspikesperbursts[:])
    spikesPerBurstsAll2.append(nspikesperbursts2[:])

    print("""forall if(ismembrane(\"sk\")) gbar_sk = gbar_sk/"""+str(gsk_apics[ig]))
    h("""forall if(ismembrane(\"sk\")) gbar_sk = gbar_sk/"""+str(gsk_apics[ig]))
    print("""forall if(ismembrane(\"cah\")) pbar_cah = pbar_cah/"""+str(gcas[ig]))
    h("""forall if(ismembrane(\"cah\")) pbar_cah = pbar_cah/"""+str(gcas[ig]))
    print("""forall if(ismembrane(\"car\")) pbar_car = pbar_car/"""+str(gcas[ig]))
    h("""forall if(ismembrane(\"car\")) pbar_car = pbar_car/"""+str(gcas[ig]))

  spTimesAllAll.append(spTimesAll[:])
  spTimesAllAll2.append(spTimesAll2[:])
  nSpikesAllAll.append(nSpikesAll[:])

#picklelist = [nSpikesAllAll,spTimesAllAll,spTimesAllAll2,condSuffixes,coeffs,Is]
#file = open('ifcurves.sav', 'w')
#pickle.dump(picklelist,file)
#file.close()
  

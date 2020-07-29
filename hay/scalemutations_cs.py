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
for idefval in range(0,len(keyList)):
  if type(defVals[keyList[idefval]]) is not list:
    defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments

unpicklefile = open('control_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
spikfreqs_control_All = unpickledlist[0]
timesc_control_All = unpickledlist[1]
Vsomac_control_All = unpickledlist[2]
VDerivc_control_All = unpickledlist[3]
VDcoeff_control_All = unpickledlist[4]
Is_control = unpickledlist[19]
Is = [0.2,0.4,0.6,0.8,1.0,1.2,1.4]

theseCoeffsAllAll = []

paramdicts = []
paramdicts.append({})                                                         # 1 spike per burst, control
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 1.6})                               # 1-2 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2})                               # 2-3 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.9})   # 3-4 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.625}) # 3-5 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.5})   # 4-6 spikes per burst      
paramdicts.append({'A_gNaTa_tbar_NaTa_t': 2.2, 'S_gCa_HVAbar_Ca_HVA': 0.3})   # 5-9 spikes per burst      

squareDurs =    [5,5,5]
targetNspikes = [[1],[1],[1]]
nextIs =        [0.0,3.0,1.5,1]
epspDivisor =   1.5
distalpoint =   620

for icell in range(0,len(paramdicts)):

  unpicklefile = open('controlamps_cs'+str(icell)+'.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  mySquareAmpsAll = unpickledlist[0]
  mySynAmpsAll = unpickledlist[1]

  squareAmps = [mySquareAmpsAll[i][0][-1] for i in range(0,len(mySquareAmpsAll))]
  synAmps = [mySynAmpsAll[i][0][-1] for i in range(0,len(mySynAmpsAll))]

  spikfreqs_control = mytools.interpolate(Is_control,spikfreqs_control_All[icell],Is)
  Vsomac_control = Vsomac_control_All[icell]
  VDerivc_control = VDerivc_control_All[icell]
  VDcoeff_control = VDcoeff_control_All[icell]
  timesc_control = timesc_control_All[icell]
  print spikfreqs_control

  morphology_file = "morphologies/cell1.asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"
  v0 = -80
  ca0 = 0.0001
  proximalpoint = 400
  distalpoint = 620
  #distalpoint = 960
  BACdt = 2.5

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
L5PC.distribute_channels("apic","gIhbar_Ih",2,-0.8696,3.6161,0.0,1.0*2.0870,0.0002)
L5PC.distribute_channels("apic","gCa_HVAbar_Ca_HVA",3,1.0,0.1,685.0,885.0,1.0*0.000555)
L5PC.distribute_channels("apic","gCa_LVAstbar_Ca_LVAst",3,1.0,0.01,685.0,885.0,1.0*0.0187)
objref sl,ns,syn1,con1,isyn, tvec
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
L5PC.apic[siteVec[0]] {
  syn1 = new AlphaSynapse(siteVec[1])
  syn1.e = 0
  syn1.tau = 5
  syn1.onset = 200 + """+str(BACdt)+""" 
  cvode.record(&syn1.i,isyn,tvec)
}
objref vsoma, vdend, vdend2, isoma, cadend, cadend2, casoma
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
L5PC.apic[siteVec[0]] cvode.record(&v(siteVec[1]),vdend2,tvec)
L5PC.apic[siteVec[0]] cvode.record(&cai(siteVec[1]),cadend2,tvec)
access L5PC.soma
isoma = new Vector()
cvode.record(&st1.i,isoma,tvec)
""")

  paramdict = paramdicts[icell]
  setparams(paramdict)

  ITERS = 20
  theseCoeffsAll = []
  theseMutValsAll = []
  theseMutVarsAll = []

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
            if mutvars[kmutvar].find('offm') > -1 or mutvars[kmutvar].find('offh') > -1 or mutvars[kmutvar].find('ehcn') > -1:
              newVal = [x+thisCoeff*mutvals for x in defVals[mutvars[kmutvar]]]
              if mutvals >= 0 and kmutvar==0:
                mutText = mutText + "+" + str(mutvals*thisCoeff) +" mV"
              elif kmutvar==0:
                mutText = mutText  + str(mutvals*thisCoeff) +" mV"
            else:
              newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvars[kmutvar]]]
              if kmutvar==0:
                mutText = mutText + "*" + str(mutvals**thisCoeff)
            if kmutvar < len(mutvars)-1:
              mutText = mutText + ", "
            #else:
            #  mutText = mutText + "\n"
            if mutvars[kmutvar].find('_Ih') > -1:
              updateThese = [1,1,1]
            elif mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_Ca_LVAst') > -1 or mutvars[kmutvar].find('_SKv3.1') > -1 or mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_SK_E2') > -1 or mutvars[kmutvar].find('_NaTa_t') > -1 or mutvars[kmutvar].find('_CaDynamics_E2') > -1:
              updateThese = [1,1,0]
            elif mutvars[kmutvar].find('_K_Pst') > -1 or mutvars[kmutvar].find('_K_Tst') > -1 or mutvars[kmutvar].find('_Nap_Et2') > -1: 
              updateThese = [1,0,0]
            elif mutvars[kmutvar].find('_Im') > -1:
              updateThese = [0,1,0]
            else:
              print "Error: str=" + str(mutvars[kmutvar])
              updatedVars = [0,0,0]
            for iupdated in range(0,3):
              if updateThese[iupdated]:
                print """forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvars[kmutvar]+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}"""
                h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvars[kmutvar]+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""")
        print mutText

        close("all")
        f, axarr = plt.subplots(2, 4)
        #for ix in range(0,3):
        #  for iy in range(0,2):
        #    axarr[iy,ix].set_position([0.05+0.3*ix, 0.05+0.4*(1-iy), 0.23, 0.3])

        tstop = 3500.0
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
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 3200
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
syn1.onset = 3200 + """+str(BACdt)+""" 
""")
            h.init()
            h.run()

            times=np.array(h.tvec)
            Vsoma=np.array(h.vsoma)
            spikes = mytools.spike_times(times,Vsoma,-20,-45)
            nSpikes1 = len(spikes)
            print "icond="+str(icond)+",iamp="+str(iampCoeff)+", nSpikes1="+str(nSpikes1)+", st1.amp="+str(h.st1.amp)+", syn1.gmax="+str(h.syn1.gmax)

            axarr[iampCoeff,icond].plot(times, Vsoma)
            axarr[iampCoeff,icond].set_title("nspikes="+str(nSpikes1))
            axarr[iampCoeff,icond].set_xlim([3190,3300])
            axarr[iampCoeff,icond].set_ylim([-100,40])

            isChanged = isChanged or nSpikes1 > 0 and iampCoeff == 0 or nSpikes1 == 0 and iampCoeff == 1


        ############################################# Condition 4: IF curve #############################################
        spikfreqs = len(Is)*[0]
        for iI in range(0,len(Is)):
          tstop = 4000.0
          squareAmp = Is[iI]
          squareDur = 3800.0
          epsp_gmax = 0.0
          h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 200
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
syn1.onset = 200 + """+str(BACdt)+""" 
""")
          h.init()
          h.run()

          times=np.array(h.tvec)
          Vsoma=np.array(h.vsoma)
          spikes = mytools.spike_times(times,Vsoma,-20,-45)
          spikfreqs[iI] = sum([1 for x in spikes if x >= 500.0])/3.5
          print "spikfreqs[iI] = "+str(sum([1 for x in spikes if x >= 500.0])/3.5)
          print "spikfreqs_control[iI] = "+str(spikfreqs_control[iI])
          if iI==4: # use the memb. pot. time course of 1.0nA for the limit cycle (not in the new version...)
            times_lc = times[:]
            Vsoma_lc = Vsoma[:]
            spikes_lc = spikes[:]

        axarr[0,3].plot(Is, spikfreqs)
        axarr[0,3].set_title("Perisomatic IF-curve")
        axarr[0,3].set_xlim([0,1.25])
        axarr[0,3].set_ylim([0,20])
        spikfreqdiffsum = sum([abs(x-y) for x,y in zip(spikfreqs,spikfreqs_control)])
        print "spikfreqdiffsum = "+str(sum([abs(x-y) for x,y in zip(spikfreqs,spikfreqs_control)]))
        spikfreqdiffrel = spikfreqdiffsum/sum(spikfreqs_control)
        print "spikfreqdiffrel = "+str(spikfreqdiffsum/sum(spikfreqs_control))

        ############################################# Condition 5: Limit cycle #############################################

        f.suptitle(mutText)
        if iter < ITERS+2:
          f.savefig("vrecs_cs"+str(icell)+"_mut"+str(igene)+"_"+str(imut)+"_"+str(iallmutval)+"_ITER"+str(iter)+".png")
        else:
          f.savefig("vrecs_cs"+str(icell)+"_mut"+str(igene)+"_"+str(imut)+"_"+str(iallmutval)+"_TEST"+str(iter-ITERS-2)+".png")

        #Print the parameters and their default values:
        for idefval in range(0,len(defVals.keys())):
          thisdefval = defVals.keys()[idefval]
          if thisdefval.find('_Im') > -1:
            h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][1]))
            #) #+" (def="+str(defVals[thisdefval])+")"
          else:
            h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][0]))
            #h('print L5PC.soma[0]."+thisdefval) #+" (def="+str(defVals[thisdefval])+")"

        isChanged = isChanged or spikfreqdiffrel > 0.1
        print isChanged
        if iter==0 and isChanged:
          print "Even null mutation causes different spiking!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
          continue
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
          if mutvars[kmutvar].find('_Ih') > -1:
            updateThese = [1,1,1]
          elif mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_Ca_LVAst') > -1 or mutvars[kmutvar].find('_SKv3.1') > -1 or mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_SK_E2') > -1 or mutvars[kmutvar].find('_NaTa_t') > -1 or mutvars[kmutvar].find('_CaDynamics_E2') > -1:
            updateThese = [1,1,0]
          elif mutvars[kmutvar].find('_K_Pst') > -1 or mutvars[kmutvar].find('_K_Tst') > -1 or mutvars[kmutvar].find('_Nap_Et2') > -1: 
            updateThese = [1,0,0]
          elif mutvars[kmutvar].find('_Im') > -1:
            updateThese = [0,1,0]
          else:
            print "Error: str=" + str(mutvars[kmutvar])
            updatedVars = [0,0,0]
          for iupdated in range(0,3):
            if updateThese[iupdated]:
              h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvars[kmutvar]+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""")
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

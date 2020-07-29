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

v0 = -80
ca0 = 0.0001
#proximalpoints = [100,100,100,200,200,200,300,300,300,400,400,400]
#distalpoints =   [600,750,900,600,750,900,600,750,900,600,750,900]
proximalpoints = [200] #,200]
distalpoints =   [600] #,900]
BACdt = 5.0
fs = 5
tstop = 5000.0
#epspdts = [0.25*x for x in range(-80,81)]
#epspdts = [0.5*x for x in range(-40,41)]
#epspdts = [2.0*x for x in range(-10,-5)]+range(-10,11)+[2.0*x for x in range(6,11)]
#epspdts = [10.0*x for x in range(-10,-8)]+[4.0*x for x in range(-20,-10)]+[2.0*x for x in range(-20,-10)]+[2.0*x for x in range(-10,11)]+[2.0*x for x in range(11,21)]+[4.0*x for x in range(11,21)]+[10.0*x for x in range(9,11)]
#epspdts_savetimecourses = [-100,-80,-60,-40,-30,-20,-10,0,10,20,30,40,60,80,100]
epspdts = [10.0*x for x in range(-10,-8)]+[4.0*x for x in range(-20,-10)]+[4.0*x for x in range(-10,-5)]+[4.0*x for x in range(-5,6)]+[4.0*x for x in range(6,11)]+[4.0*x for x in range(11,21)]+[10.0*x for x in range(9,11)]
istart = [int(100+x) for x in epspdts]
iend = [int(120+x) for x in epspdts]
#epspdts_savetimecourses = [-40,-20,-12,-4,0,4,12,20,40]
epspdts_savetimecourses = [-100,-80,-60,-32,0,32,60,80,100]

Is_st2 = 1.32
st2coeff = 0.80      #Somatic 5ms pulse
st2coeff_down = 1.35 #Somatic 5ms pulse
st1coeff = 0.9      #Proximal apical 200ms pulse

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

rateCoeffs = [[1.1],[0.7]]


counter = -1
icell = 0
idist = 0
ext2a = '_syn1coeff0.25'
ext2b = '_syn1coeff0.5'

if not exists('updownresponsemean20ms_noisy_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2a+'_controls.sav'):
 Casoma_control = []
 SKsoma_control = []
 Vdend_control = []
 Vdend2_control = []
 Vdend3_control = []
 Cadend3_control = []
 SKdend3_control = []
 spikes_control = []
 Vdend3b_control = []
 Cadend3b_control = []
 SKdend3b_control = []
 Vdend3b2_control = []
 Cadend3b2_control = []
 SKdend3b2_control = []
 for iup in range(0,2):
  if iup==0:
    ext = 'up2'
    ext2 = ext2a
  else:
    ext = 'down'
    ext2 = ext2b
            
  Casoma_control_thisUp = []
  SKsoma_control_thisUp = []
  Vdend_control_thisUp = []
  Vdend2_control_thisUp = []
  Vdend3_control_thisUp = []
  Cadend3_control_thisUp = []
  SKdend3_control_thisUp = []
  spikes_control_thisUp = []
  Vdend3b_control_thisUp = []
  Cadend3b_control_thisUp = []
  SKdend3b_control_thisUp = []
  Vdend3b2_control_thisUp = []
  Cadend3b2_control_thisUp = []
  SKdend3b2_control_thisUp = []
  for icoeff in range(0,len(rateCoeffs[iup])):
    Casoma_control_thisCoeff = []
    SKsoma_control_thisCoeff = []
    Vdend_control_thisCoeff = []
    Vdend2_control_thisCoeff = []
    Vdend3_control_thisCoeff = []
    Cadend3_control_thisCoeff = []
    SKdend3_control_thisCoeff = []
    spikes_control_thisCoeff = []
    Vdend3b_control_thisCoeff = []
    Cadend3b_control_thisCoeff = []
    SKdend3b_control_thisCoeff = []
    for rdSeed in range(1,126):
      if rdSeed%20 == 0:
        print str(100*rdSeed*0.05)+" % done of "+ext
      if not exists('updownresponse_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_0_0_0_-1_'+str(rateCoeffs[iup][icoeff])+'_seed'+str(rdSeed)+'_tmp.sav'):
        print 'updownresponse_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_0_0_0_-1_'+str(rateCoeffs[iup][icoeff])+'_seed'+str(rdSeed)+'_tmp.sav not found'
        continue
      unpicklefile = open('updownresponse_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_0_0_0_-1_'+str(rateCoeffs[iup][icoeff])+'_seed'+str(rdSeed)+'_tmp.sav', 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      spikes_control_thisCoeff.append(unpickledlist[8])
      Vdend3b_control_thisCoeff.append(unpickledlist[9])
      Cadend3b_control_thisCoeff.append(unpickledlist[10])
      SKdend3b_control_thisCoeff.append(unpickledlist[11])

    spikes_control_thisUp.append([[spikes_control_thisCoeff[i][idt] for i in range(0,len(spikes_control_thisCoeff))] for idt in range(0,len(epspdts))])
    Vdend3b_control_thisUp = [[mean(Vdend3b_control_thisCoeff[i][idt][istart[idt]:iend[idt]]) for i in range(0,len(Vdend3b_control_thisCoeff))] for idt in range(0,len(epspdts))]
    Cadend3b_control_thisUp = [[mean(Cadend3b_control_thisCoeff[i][idt][istart[idt]:iend[idt]]) for i in range(0,len(Cadend3b_control_thisCoeff))] for idt in range(0,len(epspdts))]
    SKdend3b_control_thisUp = [[mean(SKdend3b_control_thisCoeff[i][idt][istart[idt]:iend[idt]]) for i in range(0,len(SKdend3b_control_thisCoeff))] for idt in range(0,len(epspdts))]
    Vdend3b2_control_thisUp = [[mean(Vdend3b_control_thisCoeff[i][idt][100:120]) for i in range(0,len(Vdend3b_control_thisCoeff))] for idt in range(0,len(epspdts))]
    Cadend3b2_control_thisUp = [[mean(Cadend3b_control_thisCoeff[i][idt][100:120]) for i in range(0,len(Cadend3b_control_thisCoeff))] for idt in range(0,len(epspdts))]
    SKdend3b2_control_thisUp = [[mean(SKdend3b_control_thisCoeff[i][idt][100:120]) for i in range(0,len(SKdend3b_control_thisCoeff))] for idt in range(0,len(epspdts))]

  spikes_control.append(spikes_control_thisUp[:])
  Vdend3b_control.append(Vdend3b_control_thisUp[:])
  Cadend3b_control.append(Cadend3b_control_thisUp[:])
  SKdend3b_control.append(SKdend3b_control_thisUp[:])
  Vdend3b2_control.append(Vdend3b2_control_thisUp[:])
  Cadend3b2_control.append(Cadend3b2_control_thisUp[:])
  SKdend3b2_control.append(SKdend3b2_control_thisUp[:])

 picklelist = [spikes_control[:], Vdend3b_control[:], Cadend3b_control[:], SKdend3b_control[:], Vdend3b2_control[:], Cadend3b2_control[:], SKdend3b2_control[:]]
 file = open('updownresponsemean20ms_noisy_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2a+'_controls.sav','w')
 pickle.dump(picklelist,file)
 file.close()


else:
 unpicklefile = open('updownresponsemean20ms_noisy_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2a+'_controls.sav','r')
 unpickledlist = pickle.load(unpicklefile)
 unpicklefile.close()
 spikes_control = unpickledlist[0]
 Vdend3b_control = unpickledlist[1]
 Cadend3b_control = unpickledlist[2]
 SKdend3b_control = unpickledlist[3]
 Vdend3b2_control = unpickledlist[4]
 Cadend3b2_control = unpickledlist[5]
 SKdend3b2_control = unpickledlist[6]

if True:
  theseCoeffsAll = theseCoeffsAllAll[icell]
  spikesAll = []
  for idist in range(0,1):#len(proximalpoints)):
    counter = counter + 1
    proximalpoint = proximalpoints[idist]
    distalpoint = distalpoints[idist]
    fixedpoint = 700
    idist_proximal = dists.index(proximalpoint)
    idist_distal = dists.index(distalpoint)

    paramdict = paramdicts[icell]

    styles = ['b-','b-']
    ##cols = ['#00aaaa','#11cc44','#55ee00','#bbaa00','#ee6600','#ff0000', '#aa00aa','#772277','#333333']
    #cols = ['#444444','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#009999','#772277','#00cc00']
    cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
    col_control = '#2222ff'
    coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
    lw = 0.5

    mutcounter = -1

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
        mutcounter = mutcounter + 1                                                                                                                                                               
        if len(sys.argv) > 1 and int(float(sys.argv[1])) != mutcounter: # and (igene!=0 or imut!=0 or iallmutval!=0): #If 0-0-0, go a bit further anyway to load the control data
          continue

        spikesThisMutVal = []
        for icoeff in range(0,len(rateCoeffs[0])):
          for iup in range(0,2):
            if iup==0:
              ext = 'up2'
              ext2 = ext2a
            else:
              ext = 'down'
              ext2 = ext2b
            close("all")
            #f_tc,axarr_tc = subplots(2*len(epspdts_savetimecourses),4)
            #if False:
            f,axarr = subplots(len(epspdts_savetimecourses),4)
            for iter in [0,2,6,8]:
              if not exists('updownresponsemean20ms_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(iter)+'_'+str(rateCoeffs[iup][icoeff])+'.sav'):
                Casoma = []
                SKsoma = []
                Vdend = []
                Vdend2 = []
                Vdend3 = []
                Cadend3 = []
                SKdend3 = []
                spikes = []
                Vdend3b = []
                Cadend3b = []
                SKdend3b = []
                for rdSeed in range(1,126):
                  if rdSeed%20 == 0:
                    print str(100*rdSeed*0.05)+" % done of "+ext+", iter="+str(iter)
                  if not exists('updownresponse_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(iter)+'_'+str(rateCoeffs[iup][icoeff])+'_seed'+str(rdSeed)+'_tmp.sav'):
                    print 'updownresponse_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(iter)+'_'+str(rateCoeffs[iup][icoeff])+'_seed'+str(rdSeed)+'_tmp.sav not found'
                    continue
                  unpicklefile = open('updownresponse_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(iter)+'_'+str(rateCoeffs[iup][icoeff])+'_seed'+str(rdSeed)+'_tmp.sav', 'r')
                  unpickledlist = pickle.load(unpicklefile)
                  unpicklefile.close()
                
                  #picklelist = [theseCoeffsAllAll,CasomaupThisMutVal,SKsomaupThisMutVal,VdendupThisMutVal,Vdend2upThisMutVal,Vdend3upThisMutVal,Cadend3upThisMutVal,SKdend3upThisMutVal,
                  #              spikesupThisMutVal,epspdts,MT]
                  spikes.append(unpickledlist[8])
                  Vdend3b.append(unpickledlist[9])
                  Cadend3b.append(unpickledlist[10])
                  SKdend3b.append(unpickledlist[11])

                #i goes from 0 to 20 (rdSeed), idt goes from 0 to 44
                spikes_thisUp = [[spikes[i][idt] for i in range(0,len(spikes))] for idt in range(0,len(epspdts))]
                Vdend3b_thisUp = [[mean(Vdend3b[i][idt][istart[idt]:iend[idt]]) for i in range(0,len(Vdend3b))] for idt in range(0,len(epspdts))]
                Cadend3b_thisUp = [[mean(Cadend3b[i][idt][istart[idt]:iend[idt]]) for i in range(0,len(Cadend3b))] for idt in range(0,len(epspdts))]
                SKdend3b_thisUp = [[mean(SKdend3b[i][idt][istart[idt]:iend[idt]]) for i in range(0,len(SKdend3b))] for idt in range(0,len(epspdts))]
                Vdend3b2_thisUp = [[mean(Vdend3b[i][idt][100:120]) for i in range(0,len(Vdend3b))] for idt in range(0,len(epspdts))]
                Cadend3b2_thisUp = [[mean(Cadend3b[i][idt][100:120]) for i in range(0,len(Cadend3b))] for idt in range(0,len(epspdts))]
                SKdend3b2_thisUp = [[mean(SKdend3b[i][idt][100:120]) for i in range(0,len(SKdend3b))] for idt in range(0,len(epspdts))]

                picklelist = [spikes_thisUp[:], Vdend3b_thisUp[:], Cadend3b_thisUp[:], SKdend3b_thisUp[:], Vdend3b2_thisUp[:], Cadend3b2_thisUp[:], SKdend3b2_thisUp[:]]
                file = open('updownresponsemean20ms_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(iter)+'_'+str(rateCoeffs[iup][icoeff])+'.sav','w')
                pickle.dump(picklelist,file)
                file.close()
              else:
                unpicklefile = open('updownresponsemean20ms_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(iter)+'_'+str(rateCoeffs[iup][icoeff])+'.sav','r')
                unpickledlist = pickle.load(unpicklefile)
                unpicklefile.close()
                spikes_thisUp = unpickledlist[0]
                Vdend3b_thisUp = unpickledlist[1]
                Cadend3b_thisUp = unpickledlist[2]
                SKdend3b_thisUp = unpickledlist[3]
                Vdend3b2_thisUp = unpickledlist[4]
                Cadend3b2_thisUp = unpickledlist[5]
                SKdend3b2_thisUp = unpickledlist[6]
 
#              for idt in range(0,len(epspdts)):
#                spikeshist = [sum([sum([1 for x in spikes[isamp][idt] if x >= 2900+5*j and x < 2900+5*(j+1)]) for isamp in range(0,len(spikes))]) for j in range(0,60)]
#                #spikeshist = [sum([sum([1 for x in spikes[isamp][idt] if x >= 2800+5*j and x < 2800+5*(j+1)]) for isamp in range(0,len(spikes))]) for j in range(0,60)]
#                Vdend3hist = [mean([Vdend3b[i][idt][j] for i in range(0,len(spikes))]) for j in range(0,301)]
#                Cadend3hist = [mean([Cadend3b[i][idt][j] for i in range(0,len(spikes))]) for j in range(0,301)]
#                SKdend3hist = [mean([SKdend3b[i][idt][j] for i in range(0,len(spikes))]) for j in range(0,301)]
#
#                axarr[idt,0].plot([2900+5*j for j in range(0,60)], spikeshist, styles[iup], color=cols[iter])
#                axarr[idt,1].plot([2900+j for j in range(0,301)], Vdend3hist, styles[iup], color=cols[iter])
#                axarr[idt,2].plot([2900+j for j in range(0,301)], Cadend3hist, styles[iup], color=cols[iter])
#                axarr[idt,3].plot([2900+j for j in range(0,301)], SKdend3hist, styles[iup], color=cols[iter])
#
#                if iter==8:
#                  spikeshist = [sum([sum([1 for x in spikes_control[iup][icoeff][isamp][idt] if x >= 2900+5*j and x < 2900+5*(j+1)]) for isamp in range(0,len(spikes_control[iup][icoeff]))]) for j in range(0,60)]
#                  #spikeshist = [sum([sum([1 for x in spikes_control[iup][icoeff][isamp][idt] if x >= 2800+5*j and x < 2800+5*(j+1)]) for isamp in range(0,len(spikes_control[iup][icoeff]))]) for j in range(0,60)]
#                  Vdend3hist = [mean([Vdend3b_control[iup][icoeff][i][idt][j] for i in range(0,len(spikes_control[iup][icoeff]))]) for j in range(0,301)]
#                  Cadend3hist = [mean([Cadend3b_control[iup][icoeff][i][idt][j] for i in range(0,len(spikes_control[iup][icoeff]))]) for j in range(0,301)]
#                  SKdend3hist = [mean([SKdend3b_control[iup][icoeff][i][idt][j] for i in range(0,len(spikes_control[iup][icoeff]))]) for j in range(0,301)]
#
#                  axarr[idt,0].plot([2900+5*j for j in range(0,60)], spikeshist, styles[iup], color=col_control)
#                  axarr[idt,1].plot([2900+j for j in range(0,301)], Vdend3hist, styles[iup], color=col_control)
#                  axarr[idt,2].plot([2900+j for j in range(0,301)], Cadend3hist, styles[iup], color=col_control)
#                  axarr[idt,3].plot([2900+j for j in range(0,301)], SKdend3hist, styles[iup], color=col_control)
#
#              for iy in range(0,len(epspdts_savetimecourses)):
#                axarr[iy,0].set_ylim([0,16])
#                axarr[iy,0].set_ylabel(str(epspdts_savetimecourses[iy])+" ms", fontsize=7)
#                for ix in range(0,4):
#                  axarr[iy,ix].set_xlim([2900,3200])
#                  for tick in axarr[iy,ix].xaxis.get_major_ticks() + axarr[iy,ix].yaxis.get_major_ticks():
#                    tick.label.set_fontsize(fs)
#              f.savefig('updownresponsehist_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(rateCoeffs[iup][icoeff])+'.eps')
#        #    spikesThisIter.append(spikes[:])
#        #  spikesThisMutVal.append(spikesThisIter[:])
#        #spikesAll.append(spikesThisMutVal[:])



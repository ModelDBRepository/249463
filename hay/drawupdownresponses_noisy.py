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
fs = 8
tstop = 5000.0
#epspdts = [0.25*x for x in range(-80,81)]
#epspdts = [0.5*x for x in range(-40,41)]
#epspdts = [2.0*x for x in range(-10,-5)]+range(-10,11)+[2.0*x for x in range(6,11)]
#epspdts = [10.0*x for x in range(-10,-8)]+[4.0*x for x in range(-20,-10)]+[2.0*x for x in range(-20,-10)]+[2.0*x for x in range(-10,11)]+[2.0*x for x in range(11,21)]+[4.0*x for x in range(11,21)]+[10.0*x for x in range(9,11)]
#epspdts_savetimecourses = [-100,-80,-60,-40,-30,-20,-10,0,10,20,30,40,60,80,100]
epspdts = [10.0*x for x in range(-10,-8)]+[4.0*x for x in range(-20,-10)]+[4.0*x for x in range(-10,-5)]+[4.0*x for x in range(-5,6)]+[4.0*x for x in range(6,11)]+[4.0*x for x in range(11,21)]+[10.0*x for x in range(9,11)]
#epspdts_savetimecourses = [-40,-20,-12,-4,0,4,12,20,40]
epspdts_savetimecourses = [-100,-80,-60,-32,0,32,60,80,100]

Is_st2 = 1.32
st2coeff = 0.40      #Somatic 5ms pulse
st2coeff_down = 1.35 #Somatic 5ms pulse
st1coeff = 0.9      #Proximal apical 200ms pulse
syn1coeff = 0.25     #Synaptic epsp-like input

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
  exit('updownresponsemean20ms_noisy_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2a+'_controls.sav not found')
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
print 'updownresponsemean20ms_noisy_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2a+'_controls.sav loaded!'

if True:
  theseCoeffsAll = theseCoeffsAllAll[icell]
  spikesAll = []
  for idist in range(0,1):#len(proximalpoints)):
    counter = counter + 1
    proximalpoint = proximalpoints[idist]
    distalpoint = distalpoints[idist]
    fixedpoint = 700

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
          Vdend3_all = []
          Cadend3_all = []
          SKdend3_all = []
          Vdend3_2_all = []
          Cadend3_2_all = []
          SKdend3_2_all = []
          Vdend3_control = []
          for iup in range(0,2):
            Vdend3 = []
            Cadend3 = []
            SKdend3 = []
            Vdend3_2 = []
            Cadend3_2 = []
            SKdend3_2 = []
            if iup==0:
              ext = 'up2'
              ext2 = ext2a
            else:
              ext = 'down'
              ext2 = ext2b
            close("all")

            #if not exists('updownresponsemean20ms_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_0_0_0_-1_'+str(rateCoeffs[iup][icoeff])+'.sav'):
            #  print 'updownresponsemean20ms_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_0_0_0_-1_'+str(rateCoeffs[iup][icoeff])+'.sav does not exist'
            #  continue
            #unpicklefile = open('updownresponsemean20ms_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_0_0_0_-1_'+str(rateCoeffs[iup][icoeff])+'.sav','r')
            #print 'updownresponsemean20ms_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_0_0_0_-1_'+str(rateCoeffs[iup][icoeff])+'.sav loaded'
            #unpickledlist = pickle.load(unpicklefile)
            #unpicklefile.close()
            #spikes = unpickledlist[0]
            #Vdend3b = unpickledlist[1] 
            #Vdend3_control.append([mean(Vdend3b[idt]) for idt in range(0,len(epspdts))])

            for iter in [0,2,6,8]:
              Vdendmaxs_thisIter = []
              Vdendlens_thisIter = []
              if not exists('updownresponsemean20ms_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(iter)+'_'+str(rateCoeffs[iup][icoeff])+'.sav'):
                print 'updownresponsemean20ms_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(iter)+'_'+str(rateCoeffs[iup][icoeff])+'.sav does not exist'
                continue

              unpicklefile = open('updownresponsemean20ms_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(iter)+'_'+str(rateCoeffs[iup][icoeff])+'.sav','r')
              print 'updownresponsemean20ms_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(iter)+'_'+str(rateCoeffs[iup][icoeff])+'.sav loaded'
              unpickledlist = pickle.load(unpicklefile)
              unpicklefile.close()
              spikes = unpickledlist[0]
              Vdend3b = unpickledlist[1]
              Cadend3b = unpickledlist[2]
              SKdend3b = unpickledlist[3]
              Vdend3b2 = unpickledlist[4]
              Cadend3b2 = unpickledlist[5]
              SKdend3b2 = unpickledlist[6]
 
              Vdend3.append([mean(Vdend3b[idt]) for idt in range(0,len(epspdts))])
              Cadend3.append([mean(Cadend3b[idt]) for idt in range(0,len(epspdts))])
              SKdend3.append([mean(SKdend3b[idt]) for idt in range(0,len(epspdts))])
              Vdend3_2.append([mean(Vdend3b2[idt]) for idt in range(0,len(epspdts))])
              Cadend3_2.append([mean(Cadend3b2[idt]) for idt in range(0,len(epspdts))])
              SKdend3_2.append([mean(SKdend3b2[idt]) for idt in range(0,len(epspdts))])

            Vdend3_all.append(Vdend3[:])
            Cadend3_all.append(Cadend3[:])
            SKdend3_all.append(SKdend3[:])
            Vdend3_2_all.append(Vdend3_2[:])
            Cadend3_2_all.append(Cadend3_2[:])
            SKdend3_2_all.append(SKdend3_2[:])

          close("all")
          f,axarr = subplots(1,1)
          f2,axs = subplots(1,2)
          iters = [0,2,6,8]
          for iiter in range(0,len(iters)):
            iter = iters[iiter]
            axarr.plot(epspdts, Vdend3_all[0][iiter],'b-',color=cols[iter])
            axarr.plot(epspdts, Vdend3_all[1][iiter],'b--',color=cols[iter])
          #axarr.plot(epspdts, Vdend3_control[0],'b-',color=col_control)
          #axarr.plot(epspdts, Vdend3_control[1],'b--',color=col_control)
          axarr.plot(epspdts, [mean(Vdend3b_control[0][idt]) for idt in range(0,len(epspdts))],'b-',color=col_control)
          axarr.plot(epspdts, [mean(Vdend3b_control[1][idt]) for idt in range(0,len(epspdts))],'b--',color=col_control)
          axarr.set_ylim([-70,-25])
          f.savefig('updownresponsemean20ms_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.eps')

          close("all")
          f,axarr = subplots(1,1)
          axarr.set_position([0.12,0.06,0.5,0.88])
          f.text(0.04, 0.9, 'D', fontsize=27)
          iters = [0,2,6,8]
          for iiter in range(0,len(iters)):
            iter = iters[iiter]
            axarr.plot(epspdts, Vdend3_2_all[0][iiter],'b-',color=cols[iter])
            axarr.plot(epspdts, Vdend3_2_all[1][iiter],'b--',color=cols[iter])
          axarr.plot(epspdts, [mean(Vdend3b2_control[0][idt]) for idt in range(0,len(epspdts))],'b-',color=col_control)
          axarr.plot(epspdts, [mean(Vdend3b2_control[1][idt]) for idt in range(0,len(epspdts))],'b--',color=col_control)
          axarr.set_ylim([-70,-25])
          for tick in axarr.xaxis.get_major_ticks() + axarr.yaxis.get_major_ticks():
              tick.label.set_fontsize(fs)
          axarr.set_xlabel('ISI (ms)',fontsize=10)
          axarr.set_ylabel('max $V_m$ (dend) (mV)',fontsize=10)          
          f.savefig('updownresponsemean20ms_2_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+ext2+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.eps')

          Nsamp = 126
          #Nsamp = 3
          dy = 0
          axs[0].set_position([0.12,0.06,0.37,0.72])
          axs[1].set_position([0.55,0.06,0.37,0.72])
          axnew = f2.add_axes([0.17,0.85,0.75,0.09])
          for iup in range(0,2):
            if iup==0:
              ext = 'up2'
              ext2 = ext2a
              syn1coeff = 0.25
              rateCoeff = 1.1
              arrow1yplus = [0,-1,0,0,0,0,5,0,0] #red
              arrow2yplus = [0,3,5,5,0,0,0,0,0]  #black
            else:
              ext = 'down'
              ext2 = ext2b
              syn1coeff = 0.5
              rateCoeff = 0.7
              dy = -8
              arrow1yplus = [-3,-3,-3,-3,-3,0,7,0,-3]
              arrow2yplus = [0,3,12,12,0,0,0,0,0]
            for iiter in range(0,len(iters)+1):
              if iiter < len(iters):
                iter = iters[iiter]
                addition = str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'_'+str(iter)+'_'
                mycol = cols[iter]
              else:
                iter = -1
                addition = '0_0_0_-1_'
                mycol = col_control
              Vdend3tc_all = []
              for rdSeed in range(1,Nsamp):
                if not exists('updownresponsetimecourse_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_syn1coeff'+str(syn1coeff)+'_'+addition+str(rateCoeff)+'_seed'+str(rdSeed)+'.sav'):
                  print 'updownresponsetimecourse_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_syn1coeff'+str(syn1coeff)+'_'+addition+str(rateCoeff)+'_seed'+str(rdSeed)+'.sav does not exists'
                  continue
                unpicklefile = open('updownresponsetimecourse_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_syn1coeff'+str(syn1coeff)+'_'+addition+str(rateCoeff)+'_seed'+str(rdSeed)+'.sav','r')
                print 'updownresponsetimecourse_noisy'+ext+'_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_syn1coeff'+str(syn1coeff)+'_'+addition+str(rateCoeff)+'_seed'+str(rdSeed)+'.sav loaded'
                unpickledlist = pickle.load(unpicklefile)
                unpicklefile.close()
                epspdts_savetimecourses = unpickledlist[2]
                Vdend3tc = []
                for idt in range(0,len(epspdts_savetimecourses)):
                  Vdend3 = unpickledlist[1][idt][6]
                  Cadend3 = unpickledlist[1][idt][7]
                  SKdend3 = unpickledlist[1][idt][8]
                  Vdend3tc.append(Vdend3[:])
                  if idt == 7 and iup == 1 and iter == -1:
                    axnew.plot(times,Vdend3,styles[iup], color='#AAAAFF',linewidth=0.5)
                times = unpickledlist[1][idt][0]
                Vdend3tc_all.append(Vdend3tc[:])

              if len(Vdend3tc_all) > 0:
                for idt in range(0,len(epspdts_savetimecourses)):
                  axs[iup].plot(times, [mean([Vdend3tc_all[isamp][idt][it] for isamp in range(0,len(Vdend3tc_all))])+idt*50 for it in range(0,len(Vdend3tc_all[0][idt]))],styles[iup], color=mycol,linewidth=0.5)
                  if idt == 7 and iup == 1 and iter == -1:
                    axnew.plot(times, [mean([Vdend3tc_all[isamp][idt][it] for isamp in range(0,len(Vdend3tc_all))]) for it in range(0,len(Vdend3tc_all[0][idt]))],styles[iup], color=mycol,linewidth=1)
                    axnew.text(2900,-50,"ISI="+str(epspdts_savetimecourses[idt])+" ms", fontsize=8)

              #axs[iup].set_ylabel(str(epspdts[iy])+" ms", fontsize=8)                                                                                                                                                                       
              for iy in range(0,len(epspdts_savetimecourses)):
                axs[iup].text(2900,iy*50-50,"ISI="+str(epspdts_savetimecourses[iy])+" ms", fontsize=8)
                if iup == 1 and iy == 6:
                  mytools.drawarrow(axs[iup],[3000+epspdts_savetimecourses[iy]]*2, [-25+iy*50+dy+arrow1yplus[iy], -42+iy*50+dy+arrow1yplus[iy]],prc=0.8,lc='#FF0000')
                else:
                  mytools.drawarrow(axs[iup],[3000+epspdts_savetimecourses[iy]]*2, [-75+iy*50+dy+arrow1yplus[iy], -58+iy*50+dy+arrow1yplus[iy]],prc=0.8,lc='#FF0000')
                mytools.drawarrow(axs[iup],[3000]*2, [-35+iy*50+dy+arrow2yplus[iy], -52+iy*50+dy+arrow2yplus[iy]],prc=0.8,lc='#000000')

            axs[iup].set_xlim([2880,3220])
            for tick in axs[iup].xaxis.get_major_ticks() + axs[iup].yaxis.get_major_ticks():
                tick.label.set_fontsize(fs)
            axs[iup].set_ylim([-82,-66+len(epspdts_savetimecourses)*50])
            axs[iup].set_xticks([])
            axs[iup].set_yticks([])
            axs[iup].plot([3150,3150],[355,380],'k-',linewidth=2)
            axs[iup].plot([3150,3200],[355,355],'k-',linewidth=2)
            axs[iup].text(3102,367,'25 mV',fontsize=8)
            axs[iup].text(3160,360,'50 ms',fontsize=8)

          axnew.set_xlim([2880,3220])
          for tick in axnew.xaxis.get_major_ticks() + axnew.yaxis.get_major_ticks():
            tick.label.set_fontsize(fs)
          axnew.set_ylim([-75,-29])
          axnew.set_yticks([-80,-60,-40])
          axnew.set_ylabel('$V_m$ (dend)\n(mV)',fontsize=10)
          axnew.set_xticks([2900,2950,3000,3050,3100,3150,3200])
          axnew.set_xticklabels(['-100','-50','0','+50','+100','+150','+200'])
          axnew.set_xlabel('$t$ (ms)',fontsize=10)
          mytools.drawarrow(axnew,[3000+epspdts_savetimecourses[7]]*2, [-75+dy+arrow1yplus[7], -58+dy+arrow1yplus[7]],prc=0.8,lc='#FF0000')
          mytools.drawarrow(axnew,[3000]*2, [-35+dy+arrow2yplus[7], -52+dy+arrow2yplus[7]],prc=0.8,lc='#000000')
          for i in range(0,2):
            f2.text(0.04, 0.9, 'A', fontsize=27)
            f2.text(0.07, 0.74, 'B', fontsize=27)
            f2.text(0.51, 0.74, 'C', fontsize=27)

          f2.savefig('updownresponsetimecourses_noisys_cs'+str(icell)+'_dist'+str(proximalpoints[idist])+'_'+str(distalpoints[idist])+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.eps')



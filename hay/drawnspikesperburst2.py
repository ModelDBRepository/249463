from neuron import h
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
iIs = [2,4,6,8,10]
spTimesAllAll = []
spTimesAllAll2 = []
nSpikesAllAll = []
ISIs_allAll = []
fs = 8

gNaTa_apics = [1.0+x*0.1 for x in range(0,13)] + [2.2 for x in range(0,31)]
gCaHVA_somas = [1.0 for x in range(0,13)] + [0.975-0.025*x for x in range(0,31)]

condSuffixes = ['Ca_HVA','Ca_LVAst','Ih','K_Pst','K_Tst','NaTa_t','Nap_Et2','SK_E2','SKv3_1']
condNames = ['gCa_HVAbar_Ca_HVA', 'gCa_LVAstbar_Ca_LVAst', 'gIhbar_Ih', 'gImbar_Im', 'gK_Pstbar_K_Pst', 'gK_Tstbar_K_Tst', 'gNaTa_tbar_NaTa_t', 'gNap_Et2bar_Nap_Et2', 'gSK_E2bar_SK_E2', 'gSKv3_1bar_SKv3_1']

def strroundhalf(x):
  if x+0.0000001-int(x) < 0.000001:
    return str(x)
  else:
    return str(int(2*x+0.0000001)*0.5)

for icell in range(0,1):
  styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
  cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
  
  close("all")
  f,axarr = subplots(len(iIs),1)

  for iiI in range(0,len(iIs)):
    axarr[iiI].set_position([0.14,0.3+0.13*iiI,0.5,0.12])
    axarr[iiI].set_xticks([0.35+x for x in [0,6,12,22,32,42]])
    axarr[iiI].set_yticks([0,2,4,6,8])
    #axarr[iiI].grid(True, 'major', 'y', color='0.65',linestyle='-')
    axarr[iiI].set_xticklabels([])
    rect1 = matplotlib.patches.Rectangle((0.9,0), 12.05, 10, color='#ffffe4')
    rect2 = matplotlib.patches.Rectangle((12.9,0), 30.05, 10, color='#ffeeee')
    axarr[iiI].add_patch(rect1)
    axarr[iiI].add_patch(rect2)
    #if iiI != 2:
    #  axarr[iiI].set_ylabel("I = "+str(Is[iIs[iiI]]))
    axarr[iiI].text(1, 5.75, 'I = '+str(Is[iIs[iiI]])+' nA',fontsize=fs+2)
    for tick in axarr[iiI].yaxis.get_major_ticks()+axarr[iiI].xaxis.get_major_ticks():
      tick.label.set_fontsize(fs+2)
    axarr[iiI].xaxis.set_ticks_position('bottom')
    axarr[iiI].yaxis.set_ticks_position('left')


  axarr[4].text(0.2, 7.8, '*',fontsize=fs+2)
  axarr[4].text(7.2, 7.8, '*',fontsize=fs+2)
  axarr[4].text(12.2, 7.8, '*',fontsize=fs+2)
  axarr[4].text(16.2, 7.8, '*',fontsize=fs+2)
  axarr[4].text(27.2, 7.8, '*',fontsize=fs+2)
  axarr[4].text(32.2, 7.8, '*',fontsize=fs+2)
  axarr[4].text(40.2, 7.8, '*',fontsize=fs+2)

  axarr[2].set_ylabel("Number of spikes per burst\n")
  axarr[0].set_xlabel("\n\n\n\n\n\n\nHay model parameter change")
  f.savefig("test.eps")

  for ig in range(0,len(gNaTa_apics)-1):
    unpicklefile = open('spikesperburst2_'+str(ig)+'.sav', 'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    spikesPerBursts = unpickledlist[0]

    for iiI in range(0,len(iIs)):
      iI = iIs[iiI]
      I = Is[iI]
      nspikehist = []
      for nspikes in range(0,21):
        nspikehist.append(sum([1 for x in spikesPerBursts[iI][2:] if x==nspikes]))
      nspikehist = [x/sum(nspikehist) for x in nspikehist]
      maxnspikes = 0
      for inspikes in range(0,21):
        if nspikehist[inspikes] > 0:
          axarr[iiI].plot([ig,ig+0.6*nspikehist[inspikes]],[inspikes, inspikes],'b-',linewidth=2)
          maxnspikes = max(maxnspikes, inspikes)
      #axarr[iiI].set_ylim([0.5, maxnspikes+0.5])
      axarr[iiI].set_ylim([0.5, 7.5])

  for itick in range(0,3):
    axarr[0].text(0.1+6*itick, -0.5, 'gCaHVA: 100%, gNat: '+strroundhalf(100*gNaTa_apics[6*itick])+'%', {},rotation=90,fontsize=fs-2)
  for itick in range(1,4):
    axarr[0].text(0.1+12+10*itick, -0.5, 'gCaHVA: '+strroundhalf(100*gCaHVA_somas[12+10*itick])+'%, gNat: 220%', {},rotation=90,fontsize=fs-2)
  f.savefig("nSpikesPerBurstsHay.eps")


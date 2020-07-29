from neuron import h
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
from pylab import *
import mytools
import pickle
import sys
import time

v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
Is = [0.65+0.025*x for x in range(0,11)]
iIs = [4,5,6,8,10]
spTimesAllAll = []
spTimesAllAll2 = []
nSpikesAllAll = []
ISIs_allAll = []
fs = 8

gsk_apics = [1.0+x*0.25 for x in range(0,13)] + [4.0 for x in range(0,20)]
gcas = [1.0 for x in range(0,13)] + [1.05+0.05*x for x in range(0,20)]

condSuffixes = ['bk','sk','cah','car','iH','iA','kslow','na']
gNames = ['gbar','gbar','pbar','pbar','gbar','gbar','gbar','gbar']

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
    axarr[iiI].set_position([0.14,0.3+0.13*iiI,0.39,0.12])
    axarr[iiI].set_xticks([0.35+x for x in [0,6,12,18,24,30]])
    axarr[iiI].set_yticks([0,2,4,6,8])
    #axarr[iiI].grid(True, 'major', 'y', color='0.65',linestyle='-')
    axarr[iiI].set_xticklabels([])
    rect1 = matplotlib.patches.Rectangle((0.9,0), 12.05, 10, color='#ffffe4')
    rect2 = matplotlib.patches.Rectangle((12.9,0), 21.05, 10, color='#ffeeee')
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
  axarr[4].text(1.2, 7.8, '*',fontsize=fs+2)
  axarr[4].text(2.2, 7.8, '*',fontsize=fs+2)
  axarr[4].text(4.2, 7.8, '*',fontsize=fs+2)
  axarr[4].text(12.2, 7.8, '*',fontsize=fs+2)
  axarr[4].text(18.2, 7.8, '*',fontsize=fs+2)
  axarr[4].text(24.2, 7.8, '*',fontsize=fs+2)

  axarr[2].set_ylabel("Number of spikes per burst\n")
  axarr[0].set_xlabel("\n\n\n\n\n\n\nAlmog model parameter change")

  for ig in range(0,len(gsk_apics)-1):
    unpicklefile = open('spikesperburst_'+str(ig)+'.sav', 'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    spikesPerBursts = unpickledlist[0]
    spikesPerBursts2 = unpickledlist[1]
    print spikesPerBursts
    print spikesPerBursts2

    for iiI in range(0,len(iIs)):
      iI = iIs[iiI]
      I = Is[iI]
      nspikehist = []
      for nspikes in range(0,21):
        nspikehist.append(sum([1 for x in spikesPerBursts2[iI][2:] if x==nspikes]))
      nspikehist = [x/sum(nspikehist) for x in nspikehist]
      maxnspikes = 0
      for inspikes in range(0,21):
        if nspikehist[inspikes] > 0:
          axarr[iiI].plot([ig,ig+0.65*nspikehist[inspikes]],[inspikes, inspikes],'b-',linewidth=2)
          maxnspikes = max(maxnspikes, inspikes)
      #axarr[iiI].set_ylim([0.5, maxnspikes+0.5])
      axarr[iiI].set_ylim([0.5, 7.5])
      axarr[iiI].set_xlim([0,26.9])
  for itick in range(0,3):
    axarr[0].text(0.1+6*itick, -0.5, 'gca: 100%, gsk: '+strroundhalf(100*gsk_apics[6*itick])+'%', {},rotation=90,fontsize=fs-2)
  for itick in range(0,2):
    axarr[0].text(0.1+18+6*itick, -0.5, 'gca: '+strroundhalf(100*gcas[18+6*itick])+'%, gsk: 400%', {},rotation=90,fontsize=fs-2)
  f.savefig("nSpikesPerBurstsAlmog.eps")


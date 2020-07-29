import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import sys
import scipy.io
from os.path import exists

f,axarr = subplots(1,3)
for i in range(0,3):
  axarr[i].set_position([0.1+0.3*i,0.5,0.23,0.45])
  f.text(0.065+0.3*i, 0.94, chr(ord('D')+i), fontsize=22)

#Is = unique([0.34+0.0025*x for x in range(0,11)]+[0.35+0.05*x for x in range(0,22)])
Is = [0.69+0.005*x for x in range(0,43)]

styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
col_control = '#2222ff'
ispDef = 0 # Consider a local maximum above -35mV a spike only if after last spike the membrane potential came below -45mV 

icell = 0
icomb = 0
if len(sys.argv) > 1:
  icomb = int(sys.argv[1])
combs_all = [ [[0,5,1], [1,2,15], [2,4,7], [3,1,0], [5,0,0], [8,5,0], [13,2,0]],        
              [[0,5,0], [1,3,0], [2,5,1], [3,1,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],
              [[0,5,1], [1,2,15], [2,4,7], [3,1,0], [5,0,0], [8,5,0], [13,3,0]],         
              [[0,5,0], [1,3,0], [2,5,1], [3,1,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],
              [[0,5,1], [1,2,15], [2,4,7], [3,1,0], [5,0,0], [8,5,0], [13,2,0]],         
              [[0,5,0], [1,3,0], [2,5,1], [3,1,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],
              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,5,0]],         
              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,4,0]],
              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,1,0]],         
              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]],
              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,5,0]],         
              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,3,0]],
              [[0,5,1], [1,2,15], [2,4,7], [3,1,1], [5,0,0], [8,5,0], [13,0,0]],         
              [[0,5,0], [1,3,0], [2,5,1], [3,0,1], [6,3,0], [8,3,0], [12,1,1], [13,5,0]] ]
lensToStart = [150.0, 300.0, 450.0, 600.0, 650.0]
startdist = int(lensToStart[2])
currCoeff = 1.1

unpicklefile = open('ifcurvesmut2_cs'+str(icell)+'_comb'+str(icomb)+'.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
ISIsThisMutVal = unpickledlist[0]
spTimesThisMutVal = unpickledlist[1+ispDef]

#Is_control = [0.35+0.05*x for x in range(0,22)]
Is_control = [0.69+0.005*x for x in range(0,43)]

unpicklefile = open('ifcurvesmut2_cs'+str(icell)+'_0_0_0.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
spTimesThisMutVal0 = unpickledlist[1+ispDef]
nSpikes_control = [sum([1 for x in spTimesThisMutVal0[5][j] if x >= 500]) for j in range(0,len(Is_control))]

somaticIs = [-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2]
#synconductances = unique([0.000005, 0.00001, 0.000015, 0.00002, 0.000025, 0.00003, 0.000035, 0.00004, 0.000045, 0.00005, 0.000055, 0.0000025, 0.0000075, 0.0000125, 0.0000175, 0.0000225, 0.0000275, 0.0000325, 0.0000375, 0.0000425, 0.0000475, 0.0000525])
synconductances = array([0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.00008, 0.00009, 0.0001])
thrs = [[0.08, 0.18, inf], [0.08, 0.18, inf], [0.08, 0.18, inf], [0.08, 0.18, inf], [0.08, 0.18, inf], [0.08, 0.18, inf], [0.08, 0.18, inf], [0,1,2,3,inf]]


iters = [0, 2, 5, 6, 8]
for iiter in range(0,len(iters)):
  iter = iters[iiter]
  if iter==5:
    continue
  nSpikes = [sum([1 for x in spTimesThisMutVal[iiter][j] if x >= 500]) for j in range(0,len(Is))]
  axarr[0].plot(Is, [x/15.5 for x in nSpikes], styles[iter],color=cols[iter],linewidth=1)
axarr[0].plot(Is_control, [x/15.5 for x in nSpikes_control], styles[iter],color=col_control,linewidth=1)


if exists('PPIcoeffs450_cs'+str(icell)+'_0.sav'):
  print ' opening PPIcoeffs450_cs'+str(icell)+'_0.sav and PPIcoeffs_complement_450_cs'+str(icell)+'_0.sav'
  unpicklefile = open('PPIcoeffs450_cs'+str(icell)+'_0.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  PPIcoeffs_control_Almog = unpickledlist[1][4]
  print str(unpickledlist[1][4])
  ISIs_control_Almog = unpickledlist[2]
else:
  print ' PPIcoeffs450_cs'+str(icell)+'_0.sav not found'

ISIs_Almog = [10.0*x for x in range(0,51)]
iters = [0, 2, 6, 8, -1]
for iiter in range(0,len(iters)):
  iter = iters[iiter]
  if iter==5:
    continue
  unpicklefile = open('PPIcoeffs'+str(startdist)+'_cs'+str(icell)+'_comb'+str(icomb)+'_iter'+str(iter)+'.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  PPIcoeffs_Almog_thisiter = [x[:] for x in unpickledlist[1][0]]
  #ISIs_Almog = unpickledlist[2][:]
  if iter >= 0:
    axarr[1].plot(ISIs_Almog,[x[2]*currCoeff for x in PPIcoeffs_Almog_thisiter],color=cols[iter],linewidth=1.0)
axarr[1].plot(ISIs_control_Almog,[x[2]*currCoeff for x in PPIcoeffs_control_Almog],color=col_control,linewidth=1.0)

iters = [0, 2, 6, 8, -1]
coding_outputs_thismut = []
coding_outputs_control = []
for isynconductance in range(0,len(synconductances)):
  coding_outputs_thiscond = []
  coding_outputs_control_thiscond = []
  synconductance = synconductances[isynconductance]
  for iI in range(0,len(somaticIs)):
    unpicklefile = open('codings_nonprop'+str(synconductance)+'_cs'+str(icell)+'_comb'+str(icomb)+'_somaticI'+str(somaticIs[iI])+'.sav', 'r') #maybe codings/codings_nonprop'+str(synconductance)+'_cs'+str(icell)+'_comb'+str(icomb)+'.sav'
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    coding_outputs = unpickledlist[2]
    myfigs = [[],[],[],[],[]]
    myfigs_control = []
    for iplot in range(0,8):
      for iiter in range(0,len(iters)):
        iter = iters[iiter]
        if iter == -1:
          myfigs_control.append([next((i for i,x in enumerate(thrs[iplot]) if x >= coding_outputs[iiter][j][iplot])) for j in range(0,128)])
        myfigs[iiter].append([next((i for i,x in enumerate(thrs[iplot]) if x >= coding_outputs[iiter][j][iplot])) for j in range(0,128)])
    coding_outputs_thiscond.append(myfigs[:])
    if iter == -1:
      coding_outputs_control_thiscond.append(myfigs_control[:])
  coding_outputs_thismut.append(coding_outputs_thiscond[:])
  if iter == -1:
    coding_outputs_control.append(coding_outputs_control_thiscond[:])

print "Analyzing controls..."
Npatterns_control = []
for icond in range(0,len(synconductances)):
  Npatterns_control.append(mean([len(unique([sum([(4**iplot)*coding_outputs_control[icond][iI][iplot][j] for iplot in range(0,8)]) for j in range(0,128)])) for iI in range(0,len(somaticIs))]))

for iiter in range(0,len(iters)):
  iter = iters[iiter]
  if iter == 5:
    continue
  Npatterns = []
  for icond in range(0,len(synconductances)):
    Npatterns.append(mean([len(unique([sum([(4**iplot)*coding_outputs_thismut[icond][iI][iiter][iplot][j] for iplot in range(0,8)]) for j in range(0,128)])) for iI in range(0,len(somaticIs))]))
  #print str(Npatterns)                                                                                                                                                                                                                  
  print "iiter="+str(iiter)+", Npatterns="+str(mean(Npatterns))+" +- "+str(std(Npatterns))
  if iter >= 0:
    axarr[2].plot(1e6*synconductances,Npatterns,'b-',color=cols[iter])
  else:
    axarr[2].plot(1e6*synconductances,Npatterns,'k--',dashes=(1,4))
axarr[2].plot(1e6*synconductances,Npatterns_control,'b-',color=col_control)

axarr[0].set_xlabel('$I$ (nA)',fontsize=9)
#axarr[0].set_xticks([0.4,0.8,1.2])
axarr[0].set_ylabel('$f$ (Hz)',fontsize=9)
axarr[0].set_yticks([0,10,20])
axarr[0].set_xlim([0.7,0.9])
axarr[0].set_ylim([0,25])

axarr[1].set_xlabel('ISI (ms)',fontsize=9)
#axarr[1].set_xticks([0,100,200,300,400,500])
axarr[1].set_ylabel('PPI factor',fontsize=9)
axarr[1].set_yticks([0.75,1.0,1.25])
axarr[1].set_xlim([0,300])
axarr[1].set_ylim([0.6,1.4])

axarr[2].set_xlabel('Single-synapse conductance (pS)',fontsize=9)
axarr[2].set_xticks([20,40,60,80,100])
axarr[2].set_ylabel('output diversity',fontsize=9)
axarr[2].set_yticks([10,20,30])
axarr[2].set_ylim([10,39])

for i in range(0,3):
  for tick in axarr[i].yaxis.get_major_ticks()+axarr[i].xaxis.get_major_ticks():
    tick.label.set_fontsize(6)

#f.savefig("figcomb"+str(icomb)+".eps")
f.savefig("figcomb.eps")

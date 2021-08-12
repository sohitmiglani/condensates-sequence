#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize


# In[2]:


#when you have replicates, average over them
def criticalPointAveraging(tcFilename):
    betaCList=[]
    phiCList=[]
    muCList=[]
    sCList=[]
    for repNo in range(tcReps):
        thisBetaC,thisSC,thisPhiC,thisMuC=np.load(tcFilename+'_rep'+str(repNo)+'.npy')
        betaCList.append(thisBetaC)
        phiCList.append(thisPhiC)
        muCList.append(thisMuC)
        sCList.append(thisSC)

    betaC=np.mean(betaCList)
    phiC=np.mean(phiCList)
    phiStd=np.std(phiCList)
    
    tList=[1/i for i in betaCList]
    tMean=np.mean(tList)
    tStd=np.std(tList)
    
    return tMean,tStd,phiC,phiStd


# In[3]:


#load all critical data with r_+=0.5 (multiL block, L24 block, random L24)

j=0.05
tcReps=3
figureSize=[8,6]
saveBool=True

labelSize=24
pointSize=12
tickSize=18


path_l24_block='tcData/tcData_L24_block/'
path_l24_rand='tcData/tcData_L24_rand/'
path_multiL_block='tcData/tcData_multiL/'


blockSeqParams=[[8,[1,2,4]],[12,[1,2,3,6]],[16,[1,2,4]],[18,[1,3]],[20,[1,2]],[24,[1,2,3,4,6,12]]]
randNames=[str(i) for i in range(10)]

#length, domain size, Tc mean, Tc std, phiC mean, phiC std
#for random, domain size=-1
criticalData=[]
lengthSet=set()


#random sequences
for i in range(len(randNames)):
    thisSeq='L24_seq'+str(randNames[i])
    thisFilename=path_l24_rand+'tcFit_'+thisSeq+'_j'+str(j)
    tMean,tStd,phiMean,phiStd=criticalPointAveraging(thisFilename)
    criticalData.append([24,-1,tMean,tStd,phiMean,phiStd])
    
#block sequences
for theseParams in blockSeqParams:
    thisL=theseParams[0]
    lengthSet.add(thisL)
    
    for thisBlockSize in theseParams[1]:
        
        thisSeq='L'+str(thisL)+'_b'+str(thisBlockSize)
        
        if thisL==24: #average
            thisFilename=path_l24_block+'tcFit_'+thisSeq+'_j'+str(j)
            
        else:
            thisFilename=path_multiL_block+'tcFit_'+thisSeq+'_j'+str(j)
            
        tMean,tStd,phiMean,phiStd=criticalPointAveraging(thisFilename)
        criticalData.append([thisL,thisBlockSize,tMean,tStd,phiMean,phiStd])


criticalData=np.array(criticalData)


# In[4]:


#plot (Tc,PhiC) for random, block, multiL

plt.figure(figsize=(figureSize[0],figureSize[1]))

for thisL in lengthSet:
    
    thisCriticalData=criticalData[abs(criticalData[:,0]-thisL)<0.001]
    
    if thisL==24:
        randData=thisCriticalData[thisCriticalData[:,1]<0]
        blockData=thisCriticalData[thisCriticalData[:,1]>0]
        
        plt.errorbar(blockData[:,2],blockData[:,4],xerr=blockData[:,3],yerr=blockData[:,5],fmt='.',markersize=pointSize,label=str(thisL)+' (Blocks)')
        plt.errorbar(randData[:,2],randData[:,4],xerr=randData[:,3],yerr=randData[:,5],fmt='.',markersize=pointSize,label=str(thisL)+' (Scrambled)')

    else:
        plt.errorbar(thisCriticalData[:,2],thisCriticalData[:,4],xerr=thisCriticalData[:,3],yerr=thisCriticalData[:,5],fmt='.',markersize=pointSize,label=str(thisL))

plt.legend(title='Length (# of motifs)',title_fontsize=20,fontsize=20)
plt.xlabel('Critical temperature $T_\mathrm{c}$',fontsize=labelSize)
plt.ylabel('Critical density $\\phi_\mathrm{c}$',fontsize=labelSize)
plt.tick_params(axis='both',labelsize=tickSize)
plt.tight_layout()
plt.ylim(ymax=0.75)


if saveBool:
    plt.savefig('fig3_multiL_phiC.svg')


# In[9]:


#plot (Tc,PhiC) for random, block (L=24) to highlight clustering


plt.figure(figsize=(figureSize[0],figureSize[1]))

thisCriticalData=criticalData[abs(criticalData[:,0]-24)<0.001]

randData=thisCriticalData[thisCriticalData[:,1]<0]
blockData=thisCriticalData[thisCriticalData[:,1]>0]

randColor='#17becf'

plt.errorbar(randData[:,2],randData[:,4],xerr=randData[:,3],yerr=randData[:,5],fmt='.',color=randColor,markersize=10,label='Scrambled')


for i in range(len(blockData)):
    thisMarker=int(blockData[i,1])
    if i==0:
        plt.errorbar(blockData[i,2],blockData[i,4],xerr=blockData[i,3],yerr=blockData[i,5],marker='$'+str(thisMarker)+'$',color='k',markersize=12,label='Blocks')
    elif thisMarker==12:
        plt.errorbar(blockData[i,2],blockData[i,4],xerr=blockData[i,3],yerr=blockData[i,5],marker='$'+str(thisMarker)+'$',color='k',markersize=16)
    else:
        plt.errorbar(blockData[i,2],blockData[i,4],xerr=blockData[i,3],yerr=blockData[i,5],marker='$'+str(thisMarker)+'$',color='k',markersize=12)

plt.legend(title='Sequence type',title_fontsize=20,fontsize=20)
plt.xlabel('Critical temperature $T_\mathrm{c}$',fontsize=labelSize)
plt.ylabel('Critical density $\\phi_\mathrm{c}$',fontsize=labelSize)
plt.tick_params(axis='both',labelsize=tickSize)
plt.tight_layout()

if saveBool:
    plt.savefig('s_scrambledClustering.svg')


# In[8]:


#tc vs domain size


plt.figure(figsize=(figureSize[0],figureSize[1]))

thisCriticalData=criticalData[abs(criticalData[:,0]-24)<0.001]
blockData=thisCriticalData[thisCriticalData[:,1]>0]

plt.errorbar(blockData[:,1],blockData[:,2],yerr=blockData[i,3],fmt='.-',markersize=12)

plt.ylabel('Critical temperature $T_\mathrm{c}$',fontsize=labelSize)
plt.xlabel('Block size $\ell$',fontsize=labelSize)
plt.tick_params(axis='both',labelsize=tickSize)
plt.xticks([1,2,3,4,6,12])
# plt.savefig('s_domainTc.svg')
if saveBool:
    plt.savefig('s_domainTc.svg',bbox_inches = "tight")


# In[10]:


#plot tc vs L for different block sizes

#make sure same as cell 9 for fig 3
markerset=['.','_','^','s','*','o']
markersizes=[16,12,10,8,14,10]

plt.figure(figsize=(figureSize[0],figureSize[1]))


blockSizeList=[1,2,3,4,6]

for counter,thisBlockSize in enumerate(blockSizeList):
    
    thisData=criticalData[abs(criticalData[:,1]-thisBlockSize)<0.001]
    if thisBlockSize==2:
        plt.errorbar(thisData[:,0],thisData[:,2],yerr=thisData[:,3],mew=3,marker=markerset[counter],markersize=markersizes[counter],color='k',label=str(thisBlockSize))

    else:
        plt.errorbar(thisData[:,0],thisData[:,2],yerr=thisData[:,3],marker=markerset[counter],markersize=markersizes[counter],color='k',label=str(thisBlockSize))

    

ax=plt.gca()
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1],title='Block size',title_fontsize=labelSize,fontsize=labelSize)
plt.ylabel('Critical temperature $T_\mathrm{c}$',fontsize=labelSize)
plt.xlabel('Length $L$',fontsize=labelSize)
plt.tick_params(axis='both',labelsize=tickSize)
plt.tight_layout()
plt.ylim(ymin=0.8)
plt.xticks([8,10,12,14,16,18,20,22,24])


# plt.savefig('fig4_multiL_domainPlot.svg')
if saveBool:
    plt.savefig('fig4_multiL_domainPlot.svg')


# In[5]:


#stoichiometry plots

#load data
path_l24_stoichiometry='tcData/tcData_stoichiometry/'
j='0.05'

#first b3
b3Flips=['10_11_15','10_21','11_15_21','16','17_22','21','3_11_15','3_15_16_17','3_15_23','3_4_16_22','3','4_15','4','5_10_11_22','5_17','9_15_17_23']

# number of flips, mean, std
b3Summary=np.zeros([5,3])
#b3 all data: # of flips and Tc. Does not include baseline
b3Data=np.zeros([len(b3Flips),2])

#first, get baseline
b3CriticalData=criticalData[(abs(criticalData[:,0]-24)<0.001)&(abs(criticalData[:,1]-3)<0.001)]
b3Summary[0]=np.array([0,b3CriticalData[0,2],0])


#load data
for counter,thisFlip in enumerate(b3Flips):
    thisSeq='L24_b3_flip_'+thisFlip
    numFlips=thisFlip.count('_')+1
    
    
    #load reps and get mean
    tempList=[]
    for reps in range(3):
        thisData=np.load(path_l24_stoichiometry+'tcFit_'+thisSeq+'_j'+str(j)+'_rep'+str(reps)+'.npy')
        tempList.append(1/thisData[0])
    
    seqMeanTemp=np.mean(tempList)
    b3Data[counter]=np.array([numFlips,seqMeanTemp])
    
#     np.savetxt('tc_'+thisSeq+'.dat',np.array([1/thisData[0]]))

#summary statistics
for i in range(1,5):
    thisData=b3Data[b3Data[:,0]==i]
    temps=thisData[:,1]
    b3Summary[i]=np.array([i,np.mean(temps),np.std(temps)])

b3Summary[:,0]=(12+b3Summary[:,0])/24


#now load random data

# number of flips, mean, std
randSummary=np.zeros([5,3])
#rand all data: # of flips and Tc. Does not include baseline
randData=np.zeros([16,2])

#first, get baseline average over all random seqs with r_+=0.5
noFlipTcs=[]
randomSeqs=criticalData[criticalData[:,1]<0]
randSummary[0]=np.array([0,np.mean(randomSeqs[:,2]),np.std(randomSeqs[:,2])])


#load data
index=0
for pm in [[13,11],[14,10],[15,9],[16,8]]:
    for seqID in range(4):
        
        seq='L24_nP_'+str(pm[0])+'_nM_'+str(pm[1])+'_'+str(seqID)
        numFlips=pm[0]-12
        
        #load replicates and get mean Tc for each seq
        seqTempList=[]
        for reps in range(3):
            thisData=np.load(path_l24_stoichiometry+'tcFit_'+seq+'_j'+str(j)+'_rep'+str(reps)+'.npy')
            seqTempList.append(1/thisData[0])
        
        meanTemp=np.mean(seqTempList)
        randData[index]=np.array([numFlips,meanTemp])
        index+=1


# #summary statistics
for i in range(1,5):
    thisData=randData[randData[:,0]==i]
    temps=thisData[:,1]
    randSummary[i]=np.array([i,np.mean(temps),np.std(temps)])

randSummary[:,0]=(12+randSummary[:,0])/24


# In[6]:


plt.figure(figsize=(figureSize[0],figureSize[1]))
plt.errorbar(b3Summary[:,0],b3Summary[:,1],yerr=b3Summary[:,2],label='Block size = 3',color='k')
plt.errorbar(randSummary[:,0],randSummary[:,1],yerr=randSummary[:,2],label='Scrambled',color='k',fmt='--')
plt.xlabel('Motif stoichiometry $a/L$',fontsize=labelSize)
plt.ylabel('Critical temperature $T_\mathrm{c}$',fontsize=labelSize)
plt.legend(fontsize=22)
plt.xticks(np.arange(0.5,0.75,0.05))

plt.tick_params(axis='both',labelsize=tickSize)
plt.tight_layout()

if saveBool:
    plt.savefig('fig3_stoichiometry.svg')


# In[ ]:





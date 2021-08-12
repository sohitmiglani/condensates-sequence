#!/usr/bin/env python
# coding: utf-8

# 1) Load data from simulation in multicanonical ensemble.
# 2) Reweight to find P(N) in grand canonical ensemble.
# 3) Reweight in chemical potential mu to find coexistence point, where vapor and liquid peaks have equal weight.

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.optimize


# convention:
#     beta0=temp of the simulation from which the current weighting was extracted (not meaningful for uniform)
#     beta1=temp of the actual simulation being analyzed
#     betaNew=temp of the simulation we want to do next, after reweighting

# In[2]:


#begin loading data
seq='L24_b3'
linearD=30
volume=linearD**3
j=0.05
length=24


beta0=0.9267
beta1=0.9287

path='histogram/'
filename='writeRecord_multicanonical_'+seq+'_beta0_'+str(beta0)+'_beta1_'+str(beta1)+'_j'+str(j)+'.dat'
weightingFilename='preweighting_'+seq+'_beta0_'+str(beta0)+'_beta1_'+str(beta1)+'.txt'


# In[3]:


#load data


thisData=pd.read_csv(path+filename,sep=' ',header=None,names=['steps','beta','N','E'])

#find step where we reach final temp
cooledStep=0
cooledTuple=np.nonzero((thisData.iloc[:,1]==beta1).values)
if len(cooledTuple[0]>0):
    cooledStep=cooledTuple[0][0]

#thermalize
thermalizedData=thisData.iloc[cooledStep:]
thermalizedData=thermalizedData.iloc[int(len(thermalizedData)/5):]

thermalizedData['w']=np.full((thermalizedData.shape[0],1),1)
print(len(thermalizedData))

#load multicanonical weighting data. It will be just -1: -1 for uniform weighting

multicanonicalWeightingArray=np.loadtxt(path+weightingFilename,skiprows=9,delimiter=' ')
multicanonicalWeighting={}
for i in range(multicanonicalWeightingArray.shape[0]):
    multicanonicalWeighting[int(multicanonicalWeightingArray[i,0])]=multicanonicalWeightingArray[i,1]

 #load mu1 from preweighting file
weightingFileStream=open(path+weightingFilename,"r")
weightingLines=weightingFileStream.readlines()
mu1=float(weightingLines[7])

weightingFileStream.close()

#plot
fig,axs=plt.subplots(1,2,figsize=(12,4))
axs[0].hist(thermalizedData['N'],np.arange(0,1200,1))
axs[0].set_xlabel('$N$')
axs[0].set_ylabel('$P(N)$')
axs[1].plot(thermalizedData['N'])
axs[1].set_xlabel('steps')
axs[1].set_ylabel('$N$')


# In[4]:


def getReweightedData(data,beta0,mu0,beta1,mu1,nWeighting):
    
    
    dataCopy=data.copy()
    
    weightArray=np.array([nWeighting[i] for i in data['N']])
    
        
    weights=-(beta1-beta0)*data['E']+(beta1*mu1-beta0*mu0)*data['N']+np.log(weightArray)
    
    dataCopy['w']=np.exp(weights)
    
    return dataCopy
        


# In[5]:


#get data with multicanonical weighting removed
dataGCE=getReweightedData(thermalizedData,beta1,mu1,beta1,mu1,multicanonicalWeighting)


plt.figure(figsize=(8,6))
plt.hist(dataGCE['N'],weights=dataGCE['w'],bins=np.arange(0,1200,1),density=True)
plt.xlabel('$N$',fontsize=17)
plt.ylabel('$P(N)$',fontsize=17)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
# plt.savefig('fig1_pN'+str(seq)+'_beta0_'+str(beta0)+'_beta1_'+str(beta1)+'_j'+str(j)+'.png',dpi=300)


# In[6]:


#given a pdf which isn't clearly separated, return the (pdf,mu1) at a new temp which has equal weights
#define equal weights using halfway between modes as threshold


# get dilute fraction by setting threshold at halfway between 2 peaks
def overlappingDiluteFraction(data):    
    
    prelimThreshold=400
    
    binEdges=np.arange(np.min(data['N'])-0.5,np.max(data['N'])+0.5,1)
    hist,_ = np.histogram(data['N'],weights=data['w'], density=True,bins=binEdges)

    
    #find modes
    countVector=binEdges[1:]-0.5
    thresholdIndex=np.argwhere(abs(countVector-prelimThreshold)<10**-5)[0,0]

    
    #in histogram bin index
    mode1Index=np.argmax(hist[0:thresholdIndex])
    mode2Index=np.argmax(hist[thresholdIndex:])+thresholdIndex
    

    realThresholdIndex=int(round((mode1Index+mode2Index)/2))
    
    
    dilute=hist[0:realThresholdIndex]    
    diluteFraction=np.sum(dilute)
    
    return diluteFraction

    
#define a function mapping mu1 to dense fraction
def reweightedDiluteFraction(trialMu1,*paramTuple):

    #unpack the tuple containing parameters
    data,beta0,mu0,beta1=paramTuple
    
    reweightedData=getReweightedData(data,beta0,mu0,beta1,trialMu1,multicanonicalWeighting)
    reweightedFraction=overlappingDiluteFraction(reweightedData)
    
    return reweightedFraction-0.5

def getEqualWeights_overlapping(data,beta0,mu0,beta1):
    
    mu1=scipy.optimize.fsolve(reweightedDiluteFraction,mu0-0.005,args=(data,beta0,mu0,beta1))
    
    solutionData=getReweightedData(data,beta0,mu0,beta1,mu1,multicanonicalWeighting)
    return (mu1[0],solutionData)


# In[7]:


#estimate means using the midway point between modes
def getMeans_overlapping(data):
    
    binEdges=np.arange(np.min(data['N'])-0.5,np.max(data['N'])+0.5,1)
    hist,_ = np.histogram(data['N'],weights=data['w'], density=True,bins=binEdges)
    
    prelimThreshold=400
    
    #find modes
    countVector=binEdges[1:]-0.5
    thresholdIndex=np.argwhere(abs(countVector-prelimThreshold)<10**-5)[0,0]

    
    #in histogram bin index
    mode1Index=np.argmax(hist[0:thresholdIndex])
    mode2Index=np.argmax(hist[thresholdIndex:])+thresholdIndex
    

    realThresholdIndex=int(round((mode1Index+mode2Index)/2))
    
    diluteData=data[data['N']<=countVector[realThresholdIndex]]
    denseData=data[data['N']>countVector[realThresholdIndex]]
    
    diluteMean=np.average(diluteData['N'],weights=diluteData['w'])
    denseMean=np.average(denseData['N'],weights=denseData['w'])
    
    phi1=diluteMean*length/(volume)
    phi2=denseMean*length/(volume)

    
    return phi1,phi2


# In[8]:


#solve for equal weights at simulation temperature beta1
beta1MuCritical,beta1Sol=getEqualWeights_overlapping(dataGCE,beta1,mu1,beta1)

#dilute and droplet densities
phi1,phi2=getMeans_overlapping(beta1Sol)
nDilute=phi1*volume/length
nDense=phi2*volume/length

plt.figure(figsize=(11,7))

plt.hist(beta1Sol['N'],weights=beta1Sol['w'],bins=np.arange(0,1100,1),density=True)
plt.xlabel('Polymer number $N$',fontsize=22)
plt.ylabel('Probability $P(N)$',fontsize=22)
plt.tick_params(axis='both',labelsize=18)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)
ax = plt.gca()
ax.yaxis.get_offset_text().set_fontsize(18)


#dense and dilute phase
plt.axvline(nDilute,color='k',linestyle='dashed')
plt.axvline(nDense,color='k',linestyle='solid')
# plt.savefig('fig1_pN'+str(seq)+'_beta0_'+str(beta0)+'_beta1_'+str(beta1)+'_j'+str(j)+'.png',dpi=300)

#legend for dense and dilute phase
phaseY=0.008
plt.figure()
plt.plot(np.array([nDilute]),np.array([phaseY]),linestyle='dashed',color='k',label='$\\phi_{\mathrm{dilute}}=$'+str(round(phi1,3)))
plt.plot(np.array([nDense]),np.array([phaseY]),linestyle='solid',color='k',label='$\\phi_{\mathrm{dense}}=$'+str(round(phi2,3)))
plt.legend(fontsize=22)
# plt.savefig('fig1_pN_legend.png',dpi=300)

#legend for sequence label
plt.figure()
plt.title('$L=24$, $d=3$',fontsize=22 )
plt.plot(np.array([0]),np.array([0]),linestyle='solid',color='k',label='$\\phi_{\mathrm{dilute}}=$'+str(round(phi1,3)))
# plt.savefig('fig1_seqLabel.png',dpi=300)

print(beta1MuCritical)


# In[9]:


# #load data from snapshot run to establish phi threshold
# seq='L24_b3'
# linearD=30
# j=0.05
# length=24


# beta0=0.9187
# beta1=0.9187
# #phi2 at this temp: 0.779490021
# #N of 877

# path='C:/Users/bgwei/research/data/grandCanonical/figures/fig1/'
# filename='writeRecord_multicanonical_'+seq+'_beta0_'+str(beta0)+'_beta1_'+str(beta1)+'_j'+str(j)+'.dat'


# thisData=pd.read_csv(path+filename,sep=' ',header=None,names=['steps','beta','N','E'])

# #find step where we reach final temp
# cooledStep=0
# cooledTuple=np.nonzero((thisData.iloc[:,1]==beta1).values)
# if len(cooledTuple[0]>0):
#     cooledStep=cooledTuple[0][0]

# #thermalize
# thermalizedData=thisData.iloc[cooledStep:]
# thermalizedData=thermalizedData.iloc[int(len(thermalizedData)/5):]


# plt.hist(thermalizedData['N'],bins=np.arange(0,1200,1))


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize


# In[2]:


# #now fit Tc data and plot
isingBeta=0.326
isingAlpha=0.110

def getTau(phi,a,b,s,phiC,tC):
     
    
    booleanCondition=(phi<=phiC)
    tau=np.zeros(phi.shape[0])    
    for i in range(phi.shape[0]):
        
        tau[i]=scipy.optimize.fsolve(lambda tau: phiC+((-1)**booleanCondition[i])*a*tau**isingBeta+s*b*tau**(1-isingAlpha)-phi[i],0.00001)
        
    t=tC*(1-tau)
    
    return t



def fitTauFunction(data,sC,phiC,tC):
    
    popt,_=scipy.optimize.curve_fit(lambda phi,a,b:getTau(phi,a,b,sC,phiC,tC),data[:,0],data[:,1],[1,1])
    
    a,b=popt[0],popt[1]
    
    return a,b


# In[3]:


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
    muC=np.mean(muCList)
    sC=np.mean(thisSC)
    
    return betaC,phiC,muC,sC


# In[4]:


def plotSequence(filename,label,color,tcFilename=None):        

    thisData=np.load(filename+'.npy')
    thisData=thisData[thisData[:,0].argsort()]
#     thisData=thisData[1:]
    
    
    temps=1/thisData[:,0]
    
    #averaged over 3 replicates, so errorbars
    
    phi1Array=thisData[:,1]
    phi1Error=thisData[:,2]
    phi2Array=thisData[:,3]
    phi2Error=thisData[:,4]


    #phi1
    plt.errorbar(phi1Array,temps,xerr=phi1Error,marker='.',markersize=pointSize,color=color,linestyle='dashed',label=label)
    #phi2
    plt.errorbar(phi2Array,temps,xerr=phi2Error,marker='.',markersize=pointSize,color=color,linestyle='dashed')
        
     

    
    #phiC
    if tcFilename!=None:
        
        #average over replicates
        betaC,phiC,muC,sC=criticalPointAveraging(tcFilename)
        
        print(betaC)
        
        tC=1/betaC
        
        plt.plot(phiC,tC,marker='*',markersize=10,color=color)
        
        #scaling form
        allData=np.concatenate((np.array([thisData[:,1],temps]).T,np.array([thisData[:,3],temps]).T),axis=0)
        
        #only fit to bottom n data points
        n=2
        sortedTemps=np.unique(np.sort(temps))
        cutoffTemp=sortedTemps[-n]
        fittingData=allData[allData[:,1]>(cutoffTemp+0.0001)]
        lowerBound=np.min(fittingData[:,0])
        upperBound=np.max(fittingData[:,0])
        
        aFit,bFit=fitTauFunction(fittingData,sC,phiC,tC)
        scalingPhi=np.linspace(lowerBound,upperBound,20)
        scalingTau=getTau(scalingPhi,aFit,bFit,sC,phiC,tC)
        
        plt.plot(scalingPhi,scalingTau)
        

            

    


# In[5]:


#get color data
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


# In[6]:


blockSizes=[1,2,3,4,6,12]

seqList=['L24_b'+str(i) for i in blockSizes]
labelList=[str(i) for i in blockSizes]
j='0.05'
tcReps=3
phaseReps=3
saveBool=True

path_mc='binodalData/'
path_tc='tcData/'

tempLabel='Temperature $T \ (\\epsilon/k_\mathrm{B})$'
scaledTempLabel='Scaled temp. $(T-T_\mathrm{c})$ / $T_\mathrm{c}$'
scaledPhiLabel='Scaled density $(\\phi-\\phi_\mathrm{c})$ / $\\phi_\mathrm{c}$'
labelSize=24
tickSize=18
figSize=[8,6]
pointSize=12

fig = plt.figure(figsize=(figSize[0]+4,figSize[1]))
ax = fig.add_subplot(1, 1, 1)

for counter,seq in enumerate(seqList):
    
    thisFilename=path_mc+'phaseData_'+seq+'_j'+str(j)
    
    thisLabel=labelList[counter]
    
    print(seq)
    
    thisTcFilename=path_tc+'tcFit_'+seq+'_j'+str(j)
    plotSequence(thisFilename,thisLabel,colors[counter],thisTcFilename)

    
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], title='Block size',title_fontsize=labelSize,fontsize=22)
plt.xlim(xmax=1)
plt.xticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])

plt.xlabel('Density $\\phi$',fontsize=labelSize)
plt.ylabel(tempLabel,fontsize=labelSize)
plt.tick_params(axis='both',labelsize=tickSize)
# plt.title('Simulations',fontsize=20)
if saveBool==True:
    plt.savefig('fig2_phaseDiagram_simulation.svg')


# In[7]:


1.3103921542214765/1.0888366012595907


# In[8]:


#data collapse hypothesis: what if we rescale by tc, phiC?
def plotSequenceScaled(filename,label,color,tcFilename):
    
    #load raw data
    thisData=np.load(filename)
    thisData=thisData[thisData[:,0].argsort()]
#     thisData=thisData[1:]
    temps=1/thisData[:,0]

    
    #load critical data
    betaC,phiC,muC,sC=criticalPointAveraging(tcFilename)
    tC=1/betaC
    
    rescaledTemp=(temps-tC)/tC
    
    rescaledPhi1=(thisData[:,1]-phiC)/phiC    
    rescaledPhi2=(thisData[:,3]-phiC)/phiC
    
    relativeError_phi1=np.divide(thisData[:,2],thisData[:,1])
    relativeError_phi2=np.divide(thisData[:,4],thisData[:,3])
    
#     scaledError_phi1=np.multiply(rescaledPhi1,relativeError_phi1)
#     scaledError_phi2=np.multiply(rescaledPhi2,relativeError_phi2)

    scaledError_phi1=thisData[:,2]/phiC
    scaledError_phi2=thisData[:,4]/phiC
    
    
    #phi1
    plt.errorbar(rescaledPhi1,rescaledTemp,xerr=scaledError_phi1,marker='.',markersize=pointSize,color=color,linestyle='dashed',label=label)
    
    #phi2
    plt.errorbar(rescaledPhi2,rescaledTemp,xerr=scaledError_phi2,marker='.',markersize=pointSize,color=color,linestyle='dashed')
    
    #critical point
    plt.plot(0,0,marker='*',markersize=10,color='k')
            

    


# In[9]:


#collapsed data

plt.figure(figsize=(figSize[0]+1,figSize[1]))
ax=plt.gca()

for counter,seq in enumerate(seqList):
    
    thisFilename=path_mc+'phaseData_'+seq+'_j'+str(j)+'.npy'
    

    thisLabel=labelList[counter]
    thisTcFilename=path_tc+'tcFit_'+seq+'_j'+str(j)
    plotSequenceScaled(thisFilename,thisLabel,colors[counter],thisTcFilename)

      
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles[::-1], labels[::-1], title='Domain size',title_fontsize=19,fontsize=19)

plt.xlabel(scaledPhiLabel,fontsize=labelSize)
plt.ylabel(scaledTempLabel,fontsize=labelSize)
plt.tick_params(axis='both',labelsize=tickSize)
# plt.ticklabel_format(style='sci',axis='y', scilimits=(0,0))
plt.ticklabel_format(axis='y', scilimits=(0,0),useMathText=True)

ax = plt.gca()
ax.yaxis.get_offset_text().set_fontsize(tickSize+5)
if saveBool==True:
    plt.savefig('fig2_phaseDiagram_simulation_scaled.svg')


# In[13]:


#plot binodal from theory

def readBinodalData(filename):
    
    with open(filename,"r") as f:
        lines=f.readlines()
    
    data=[]
    for thisLine in lines:
        if thisLine.find('none')==-1:
            lineArray=np.fromstring(thisLine,sep=' ')
            data.append(lineArray)
        
    binodalArray=np.array(data)
    return binodalArray

#colors
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

seqList=['L24_b'+str(i) for i in [1,2,3,4,6,12]]
labelList=[str(i) for i in [1,2,3,4,6,12]]


path_theory='C:/Users/bgwei/research/data/grandCanonical/figures/fig2/theory/'

fig = plt.figure(figsize=(figSize[0]+2,figSize[1]))
ax = fig.add_subplot(1, 1, 1)

for counter,seq in enumerate(seqList):
    
    tcData=np.loadtxt(path_theory+'criticalPoint_'+seq+'_j-'+j+'.dat')
    print(seq,1/tcData[1])
    
    binodalFilename=path_theory+'theoryBinodal_'+seq+'_j'+j+'.dat'
    binodalData=readBinodalData(binodalFilename)
    
        
    thisLabel=labelList[counter]
    
    #add tc to binodalData to connect lines
    tcPoint=np.array([[0,tcData[1],tcData[0],tcData[0]]])
    binodalData=np.concatenate((tcPoint,binodalData),axis=0)
    
    
    temps=1/binodalData[:,1]
    
    #reveal blue curve
    if seq=='L24_b2':
        binodalData=binodalData[::2]
        temps=temps[::2]
    
    #phi1
    plt.plot(binodalData[:,2],temps,marker='.',markersize=pointSize,color=colors[counter],linestyle='dashed',label=thisLabel)
    
    #phi2
    plt.plot(binodalData[:,3],temps,marker='.',markersize=pointSize,color=colors[counter],linestyle='dashed')
    
    #cp star
    plt.plot(tcData[0],1/tcData[1],marker='*',markersize=10,color=colors[counter])


# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles[::-1], labels[::-1], title='Domain size',title_fontsize=15,fontsize=15)
plt.xlabel('Density $\\phi$',fontsize=labelSize)
plt.ylabel(tempLabel,fontsize=labelSize)
# plt.title('Mean field theory',fontsize=20)
plt.tick_params(axis='both',labelsize=tickSize)
plt.title('Mean-field theory',fontsize=labelSize)
# plt.savefig('fig2_phaseDiagram_theory.svg')
if saveBool==True:
    plt.savefig('fig2_phaseDiagram_theory.svg')


# In[11]:


1.4439699912252202/1.2991368501958858


# In[14]:


#plot binodal from theory, rescaled

def readScaledBinodalData(filename,phiC):
    
    with open(filename,"r") as f:
        lines=f.readlines()
    
    data=[]
    for thisLine in lines:
        if thisLine.find('none')==-1:
            lineArray=np.fromstring(thisLine,sep=' ')
            
            lineArray[2]=(lineArray[2]-phiC)/phiC
            lineArray[3]=(lineArray[3]-phiC)/phiC

            data.append(lineArray)
        
    binodalArray=np.array(data)
    return binodalArray

#colors
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

seqList=['L24_b'+str(i) for i in [1,2,3,4,6,12]]

j='0.05'

path_theory='C:/Users/bgwei/research/data/grandCanonical/figures/fig2/theory/'


plt.figure(figsize=(figSize[0]+2,figSize[1]))

for counter,seq in enumerate(seqList):
    
    tcData=np.loadtxt(path_theory+'criticalPoint_'+seq+'_j-'+j+'.dat')
    
    binodalFilename=path_theory+'theoryBinodal_'+seq+'_j'+j+'.dat'
    binodalData=readScaledBinodalData(binodalFilename,tcData[0])
        
    #add tc to binodalData to connect lines
    tcPoint=np.array([[0,tcData[1],0,0]])
    binodalData=np.concatenate((tcPoint,binodalData),axis=0)    
        
    
    t=np.power(binodalData[:,1],-1)
    tau=(t-1/tcData[1])/(1/tcData[1])
    
    #phi1
    plt.plot(binodalData[:,2],tau,marker='.',markersize=10,color=colors[counter],linestyle='dashed',label=blockSizes[counter])
    
    #phi2
    plt.plot(binodalData[:,3],tau,marker='.',markersize=10,color=colors[counter],linestyle='dashed')
    
    plt.plot(np.array([0]),np.array([0]),marker='*',markersize=10,color='k')

plt.xlabel(scaledPhiLabel,fontsize=labelSize)
plt.ylabel(scaledTempLabel,fontsize=labelSize)
plt.tick_params(axis='both',labelsize=tickSize)
plt.title('Mean-field theory',fontsize=labelSize)

# plt.legend(title='Domain size',title_fontsize=15,fontsize=15)

# plt.savefig('fig2_phaseDiagram_theory_scaled.png',dpi=300)
# plt.savefig('fig2_phaseDiagram_theory_scaled.svg')

if saveBool==True:
    plt.savefig('fig2_phaseDiagram_theory_scaled.svg')
    


# In[ ]:




